"""Verify PR 22's UserWarning path survives PR 23's strict pytest config.

PR 23 adds filterwarnings = error::RuntimeWarning. PR 22's
resolve_r_form() emits warnings.warn(..., UserWarning) on the
cw_validated=False path (plywood). The narrow RuntimeWarning filter
must NOT convert UserWarning into errors; this test exercises the
plywood path and asserts the UserWarning is emitted (not raised).
"""
import warnings

import pytest

from thermal_engine_2d import _warned_form_types, resolve_r_form


class _Stub:
    def __init__(self, form_type):
        self.form_type = form_type


def test_resolve_r_form_steel_no_warning():
    _warned_form_types.clear()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = resolve_r_form(_Stub("steel"))
    assert r == pytest.approx(0.0862)
    assert len(w) == 0


def test_resolve_r_form_plywood_emits_userwarning():
    _warned_form_types.clear()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = resolve_r_form(_Stub("plywood"))
    assert r == pytest.approx(0.17)
    assert len(w) == 1
    assert issubclass(w[0].category, UserWarning)
    assert "out-of-envelope" in str(w[0].message).lower()


def test_resolve_r_form_plywood_warning_emitted_once():
    _warned_form_types.clear()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        resolve_r_form(_Stub("plywood"))
        resolve_r_form(_Stub("plywood"))
    assert len(w) == 1


def test_resolve_r_form_plastic_liner_raises():
    _warned_form_types.clear()
    with pytest.raises(NotImplementedError):
        resolve_r_form(_Stub("plastic_liner"))


def test_resolve_r_form_unknown_raises_keyerror():
    _warned_form_types.clear()
    with pytest.raises(KeyError) as exc_info:
        resolve_r_form(_Stub("titanium"))
    msg = str(exc_info.value)
    assert "steel" in msg
    assert "plywood" in msg
