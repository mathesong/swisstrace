from datetime import datetime, timedelta, timezone
import numpy as np
import pytest


@pytest.fixture
def raw_factory(tmp_path):
    def write(counts=None, times=None, filename="raw.crv", start=None):
        if times is None:
            times = np.arange(301, dtype=float)
        if counts is None:
            counts = np.full(len(times), 30.)
            post = times >= 40
            counts[post] += 100 * np.exp(-np.log(2) / 6586.2 * (times[post] - 40))
        start = start or datetime(2026, 1, 1, tzinfo=timezone.utc)
        path = tmp_path / filename
        with path.open("w") as stream:
            for t, count in zip(times, counts):
                stamp = start + timedelta(seconds=float(t))
                stream.write(f'{stamp:%Y %m %d %H %M} {stamp.second + stamp.microsecond / 1e6:.6f} {count:.17g} 550 1000\n')
        return path
    return write
