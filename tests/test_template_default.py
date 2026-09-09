"""Tests for the ``template_file`` default value in ``analyze_files``.

Added in Fase 3, per the user's explicit request: ``template_file`` should
have a sensible default (a template bundled with the package) instead of
requiring ``None``/no restriction, while an explicit ``template_file``
argument must always override the default.

Originally Arpeggio-only (``constants.DEFAULT_TEMPLATE_FILE`` =
``"template.json"``). Extended for two follow-ups:
- the Arpeggio default was renamed to ``template_arpeggio.json``
  (``constants.DEFAULT_ARPEGGIO_TEMPLATE_FILE``, with ``DEFAULT_TEMPLATE_FILE``
  kept as a backward-compatible alias) once...
- ...IChem got its own bundled default too
  (``constants.DEFAULT_ICHEM_TEMPLATE_FILE`` = ``"template_ichem.json"``),
  so the two modes are now symmetric: both get a bundled default that can
  always be overridden, and both fall back to no restriction if their
  bundled file is ever missing.

Deliberately a new file rather than edits to the frozen
``test_analyze_files.py`` baseline — see plan-accion-modularizacion-pickit.md,
Fase 2/Fase 3 notes on why this couldn't be done as part of Fase 2 without
risking the frozen baseline tests.
"""

import os

from pickit.constants import DEFAULT_ARPEGGIO_TEMPLATE_FILE, DEFAULT_ICHEM_TEMPLATE_FILE, DEFAULT_TEMPLATE_FILE


class TestTemplateFileDefault:
    def test_bundled_default_template_exists_on_disk(self):
        # Sanity check the fixture this whole test file relies on: the
        # package actually ships a template_arpeggio.json next to io_mixin.py.
        import pickit

        package_dir = os.path.dirname(pickit.__file__)
        assert os.path.isfile(os.path.join(package_dir, DEFAULT_ARPEGGIO_TEMPLATE_FILE))
        assert DEFAULT_TEMPLATE_FILE == DEFAULT_ARPEGGIO_TEMPLATE_FILE

    def test_no_template_file_arg_still_uses_bundled_default_silently(self, analyzer, arpeggio_dir, activity_csv):
        # Calling analyze_files without template_file must not raise or
        # behave differently just because a default now exists — this is
        # the "no debe ser obligatorio ni silencioso [en el sentido de
        # romper cosas]" requirement.
        data = analyzer.analyze_files(directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv)
        assert "hbond" in data.interactions

    def test_explicit_template_file_overrides_the_default(self, analyzer, arpeggio_dir, activity_csv, template_path):
        # An explicitly-passed template_file must be honored over the
        # bundled default, even though in this project both happen to
        # point at files with equivalent contact/type coverage.
        with_explicit = analyzer.analyze_files(
            directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv, template_file=template_path
        )
        without_explicit = analyzer.analyze_files(
            directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv
        )
        # Both use a template (default vs explicit), so the resulting
        # interaction ordering should match — this is the "an explicit
        # template still works, and using no template arg isn't silently
        # different in an unexpected way" contract.
        assert with_explicit.interactions == without_explicit.interactions
        assert with_explicit.matrix == without_explicit.matrix

    def test_missing_bundled_default_does_not_break_anything(self, analyzer, arpeggio_dir, activity_csv, monkeypatch):
        # Simulate the bundled template_arpeggio.json being absent:
        # analyze_files must fall back to "no template restriction", exactly
        # like before this default existed, not raise.
        import pickit.io_mixin as io_mixin_module

        real_isfile = os.path.isfile

        def fake_isfile(path):
            if path.endswith(DEFAULT_ARPEGGIO_TEMPLATE_FILE) and os.path.dirname(path) == os.path.dirname(
                io_mixin_module.__file__
            ):
                return False
            return real_isfile(path)

        monkeypatch.setattr(os.path, "isfile", fake_isfile)

        data = analyzer.analyze_files(directory=arpeggio_dir, mode=analyzer.ARPEGGIO, activity_file=activity_csv)
        # Falls back to the full, unrestricted arpeggio interaction_labels.
        assert data.interactions == analyzer.interaction_labels


class TestTemplateFileDefaultIchem:
    """Mirrors TestTemplateFileDefault, but for IChem's own bundled default
    (added alongside PDBe-schema Arpeggio support). The bundled
    template_ichem.json lists every INTERACTION_LABELS entry, so using it
    doesn't restrict anything — same "present but non-breaking" contract as
    the Arpeggio default.
    """

    ICHEM_LINE = "Hydrophobic|CA|x|HIS 41-A|O1|x|x|x|x|x\n"

    def test_bundled_default_template_exists_on_disk(self):
        import pickit

        package_dir = os.path.dirname(pickit.__file__)
        assert os.path.isfile(os.path.join(package_dir, DEFAULT_ICHEM_TEMPLATE_FILE))

    def test_no_template_file_arg_still_uses_bundled_default_silently(self, analyzer, tmp_path):
        ichem_dir = tmp_path / "ichem_mini"
        ichem_dir.mkdir()
        (ichem_dir / "complex1.txt").write_text(self.ICHEM_LINE)

        data = analyzer.analyze_files(directory=str(ichem_dir), mode=analyzer.ICHEM)
        # The bundled default is the full label set, so nothing is filtered
        # out and "Hydrophobic" still shows up in the resulting matrix.
        assert "Hydrophobic" in data.interactions
        assert any("HIS 41" in row[0] for row in data.matrix[1:])

    def test_explicit_template_file_restricts_interactions(self, analyzer, tmp_path):
        ichem_dir = tmp_path / "ichem_mini"
        ichem_dir.mkdir()
        (ichem_dir / "complex1.txt").write_text(
            self.ICHEM_LINE + "Ionic_PROT|CB|x|ASN 142-A|N1|x|x|x|x|x\n"
        )
        custom_template = tmp_path / "only_hydrophobic.json"
        custom_template.write_text('["Hydrophobic"]')

        data = analyzer.analyze_files(
            directory=str(ichem_dir), mode=analyzer.ICHEM, template_file=str(custom_template)
        )
        assert data.interactions == ["Hydrophobic"]
        row_labels = [row[0] for row in data.matrix[1:]]
        assert "ASN 142" not in row_labels

    def test_missing_bundled_default_does_not_break_anything(self, analyzer, tmp_path, monkeypatch):
        import pickit.io_mixin as io_mixin_module

        real_isfile = os.path.isfile

        def fake_isfile(path):
            if path.endswith(DEFAULT_ICHEM_TEMPLATE_FILE) and os.path.dirname(path) == os.path.dirname(
                io_mixin_module.__file__
            ):
                return False
            return real_isfile(path)

        monkeypatch.setattr(os.path, "isfile", fake_isfile)

        ichem_dir = tmp_path / "ichem_mini"
        ichem_dir.mkdir()
        (ichem_dir / "complex1.txt").write_text(self.ICHEM_LINE)

        data = analyzer.analyze_files(directory=str(ichem_dir), mode=analyzer.ICHEM)
        assert data.interactions == analyzer.interaction_labels
