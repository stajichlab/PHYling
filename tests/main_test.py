from __future__ import annotations

import pytest

import phyling.__main__ as main
from phyling import VERSION


def test_main_help():
    assert main.main(["--help"]) == main.main(["-h"]) == 0


def test_main_version(capsys: pytest.CaptureFixture):
    assert main.main(["--version"]) == main.main(["-V"])
    captured: str = capsys.readouterr().out
    captured = captured.split("\n")
    assert captured[0] == captured[1] == VERSION
