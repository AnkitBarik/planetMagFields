import logging

import pytest

from planetmagfields.cli import main


@pytest.mark.parametrize("argv", [
    ["-p", "jupiter", "-o", "nosuchmodel"],   # missing data file
    ["-p", "pluto"],                          # unknown planet
])
def test_error_reported_once(argv, caplog):
    with caplog.at_level(logging.INFO):
        assert main(argv) == 1
    errors = [r for r in caplog.records if r.levelno >= logging.ERROR]
    assert len(errors) == 1, [r.getMessage() for r in errors]
