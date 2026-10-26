import pandas as pd
import pytest

from plate_truth.report import (
    fmt_count,
    fmt_interval,
    fmt_pct,
    fmt_share,
    group_habitats,
    level_presence,
    omitted_level,
    study_label,
)


def test_study_label():
    assert study_label("indraningrat2019") == "Indraningrat 2019"
    with pytest.raises(ValueError):
        study_label("zheng2024_soil")


def test_group_habitats_joins_in_order_without_repeats():
    census = pd.DataFrame({"study_group": ["z", "z", "z", "k"],
                           "habitat": ["soil", "sludge", "soil", "soil"]})
    assert group_habitats(census) == {"z": "soil + sludge", "k": "soil"}


def test_level_presence():
    ls = pd.DataFrame({"factor": ["medium", "oxygen"], "level": ["rich", "anaerobic"],
                       "studies": ["a;b", "b"]})
    out = level_presence(ls, ["a", "b", "c"])
    assert out.loc["a"].tolist() == [True, False]
    assert out.loc["b"].tolist() == [True, True]
    assert out.loc["c"].tolist() == [False, False]
    assert list(out.columns) == ["medium:rich", "oxygen:anaerobic"]


def test_omitted_level():
    est = pd.Series({"m_a": 0.3, "m_b": -0.1, "x": 5.0})
    v = pd.DataFrame([[0.04, 0.01, 0.0], [0.01, 0.09, 0.0], [0.0, 0.0, 1.0]],
                     index=est.index, columns=est.index)
    e, se = omitted_level(est, v, ["m_a", "m_b"])
    assert e == pytest.approx(-0.2)
    assert se == pytest.approx((0.04 + 0.09 + 2 * 0.01) ** 0.5)


def test_formatting():
    assert fmt_share(0.5012) == "0.50"
    assert fmt_pct(0.0331) == "3.3%"
    assert fmt_interval(0.35, 0.59) == "(0.35-0.59)"
    assert fmt_interval(-0.07, 0.36) == "(-0.07 to 0.36)"
    assert fmt_interval(0.011, 0.22, digits=1, pct=True) == "(1.1%-22.0%)"
    assert fmt_count(665023) == "665,023"
    with pytest.raises(ValueError):
        fmt_count(2.5)
