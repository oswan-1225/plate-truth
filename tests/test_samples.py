import pandas as pd

from plate_truth.samples import condition_id, parse_zheng_library


def test_parse_zheng_library_kinds():
    names = pd.Series(["TA_3_ORI_ORI", "TA_7_FIL_FIL", "TA_12_AN_TSA", "GD_1_MI_inSITU",
                       "XM2_AE_31_FIL_FIL"])
    out = parse_zheng_library(names)
    assert out["kind"].tolist()[:4] == ["environment", "inoculum", "plate", "plate"]
    assert out.loc[2, "medium"] == "1% TSA" and out.loc[2, "oxygen"] == "anaerobic"
    assert out.loc[3, "medium"] == "soil extract" and out.loc[3, "oxygen"] == "microaerophilic"
    assert pd.isna(out.loc[0, "medium"])  # environment rows carry no medium
    assert out.loc[2, "replicate"] == 12
    assert out.iloc[4].isna().all()  # sludge-style name does not fit the soil pattern


def test_condition_id_only_for_plates():
    df = parse_zheng_library(pd.Series(["SM_1_AE_TSA", "SM_2_ORI_ORI"]))
    cid = condition_id(df)
    assert cid[0] == "SM|1% TSA|aerobic"
    assert pd.isna(cid[1])
