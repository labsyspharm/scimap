import importlib
from unittest.mock import Mock

import matplotlib
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

matplotlib.use("Agg")


def _plotting_adata():
    obs = pd.DataFrame(
        {
            "imageid": ["one", "one", "two", "two", "two", "two"],
            "phenotype": pd.Categorical(["A", "A", "B", "B", "C", "C"]),
            "X_centroid": [0, 2, 0, 2, 1, 3],
            "Y_centroid": [0, 0, 2, 2, 1, 1],
            "zone": ["inner", "outer", "inner", "outer", "inner", "outer"],
        },
        index=[f"cell-{i}" for i in range(6)],
    )
    adata = AnnData(
        np.arange(18, dtype=float).reshape(6, 3),
        obs=obs,
        var=pd.DataFrame(index=["M1", "M2", "M3"]),
    )
    adata.raw = adata.copy()
    adata.layers["scaled"] = adata.X / 10
    adata.obsm["umap"] = np.arange(12, dtype=float).reshape(6, 2)
    return adata


def _interaction_data():
    rows = []
    for phenotype in ["A", "B"]:
        for neighbour in ["A", "B"]:
            rows.append(
                {
                    "phenotype": phenotype,
                    "neighbour_phenotype": neighbour,
                    "one": 2 if phenotype == neighbour else -1,
                    "two": 4 if phenotype == neighbour else -2,
                    "pvalue_one": 0.01 if phenotype == neighbour else 0.2,
                    "pvalue_two": 0.02 if phenotype == neighbour else 0.3,
                }
            )
    return pd.DataFrame(rows)


def _distance_data(adata):
    return pd.DataFrame(
        {
            "A": [0.5, 1.0, 2.0, 2.5, 3.0, 3.5],
            "B": [1.5, 1.0, 0.5, 1.0, 1.5, 2.0],
            "C": [3.0, 2.5, 2.0, 1.5, 1.0, 0.5],
        },
        index=adata.obs_names,
    )


def test_umap_returns_uncoloured_and_mixed_colour_data(monkeypatch, tmp_path):
    module = importlib.import_module("scimap.plotting.umap")

    adata = _plotting_adata()
    monkeypatch.setattr(module.plt, "show", Mock())

    plain = module.umap(adata, return_data=True, tight_layout=True)
    assert plain.columns.tolist() == ["umap-1", "umap-2"]

    mixed = module.umap(
        adata,
        color=["phenotype", "M1"],
        layer="scaled",
        palette={"A": "black", "B": "red", "C": "blue"},
        ncols=2,
        return_data=True,
        saveDir=tmp_path,
        fileName="mixed.png",
    )
    assert mixed.columns.tolist() == ["umap-1", "umap-2", "phenotype", "M1"]
    np.testing.assert_allclose(mixed["M1"], adata.layers["scaled"][:, 0])
    assert (tmp_path / "mixed.png").is_file()


def test_umap_uses_raw_and_validates_requested_data(monkeypatch):
    module = importlib.import_module("scimap.plotting.umap")

    adata = _plotting_adata()
    monkeypatch.setattr(module.plt, "show", Mock())

    returned = module.umap(adata, color="M3", use_raw=True, return_data=True)
    np.testing.assert_allclose(returned["M3"], adata.raw.X[:, 2])

    with pytest.raises(ValueError, match="not found"):
        module.umap(adata, color="not-a-marker")
    del adata.obsm["umap"]
    with pytest.raises(KeyError, match="run `sm.tl.umap"):
        module.umap(adata)


def test_heatmap_exercises_manual_and_clustered_modes(monkeypatch, tmp_path):
    module = importlib.import_module("scimap.plotting.heatmap")

    adata = _plotting_adata()
    monkeypatch.setattr(module.plt, "show", Mock())

    module.heatmap(
        adata,
        groupBy="phenotype",
        layer="raw",
        subsetMarkers=["M3", "M1"],
        clusterRows=False,
        clusterColumns=False,
        standardScale="column",
        orderRow=["C", "B", "A"],
        orderColumn=["M1", "M3"],
        showPrevalence=True,
        saveDir=tmp_path,
        fileName="ordered.png",
        vmin=-2,
        vmax=2,
    )
    assert (tmp_path / "ordered.png").is_file()

    module.heatmap(
        adata,
        groupBy="phenotype",
        layer="scaled",
        clusterRows=True,
        clusterColumns=True,
        standardScale="row",
    )
    with pytest.raises(ValueError, match="Both 'saveDir' and 'fileName'"):
        module.heatmap(adata, groupBy="phenotype", fileName="orphan.png")
    with pytest.raises(ValueError, match="standardScale"):
        module.heatmap(
            adata,
            groupBy="phenotype",
            clusterRows=False,
            clusterColumns=False,
            standardScale="invalid",
        )


def test_voronoi_supports_limits_overlay_and_custom_colours(monkeypatch, tmp_path):
    module = importlib.import_module("scimap.plotting.voronoi")

    adata = _plotting_adata()
    monkeypatch.setattr(module.plt, "show", Mock())

    module.voronoi(
        adata,
        subset="two",
        color_by="phenotype",
        colors={"A": "red", "B": "green", "C": "blue"},
        x_lim=[0, 3],
        y_lim=[0, 2],
        flip_y=False,
        voronoi_edge_color="facecolor",
        voronoi_alpha=[0.2, 0.4, 0.6, 0.8],
        overlay_points="zone",
        overlay_points_categories="inner",
        overlay_points_colors={"inner": "purple"},
        overlay_point_shape="x",
        saveDir=tmp_path,
        fileName="voronoi.png",
    )
    assert (tmp_path / "voronoi.png").is_file()

    with pytest.raises(ValueError, match="Color mapping"):
        module.voronoi(
            adata,
            subset="two",
            color_by="phenotype",
            colors={"B": "green"},
            plot_legend=False,
        )


def test_add_roi_omero_combines_overlaps_and_assigns_buffer_borders(capsys):
    from scimap.helpers.addROI_omero import addROI_omero

    obs = pd.DataFrame(
        {
            "X_centroid": [3.0, 7.0, -1.0, 13.0, 50.0],
            "Y_centroid": [5.0, 5.0, 5.0, 5.0, 5.0],
            "imageid": ["one", "one", "one", "one", "two"],
        },
        index=["alpha", "overlap", "border", "beta", "other-image"],
    )
    adata = AnnData(np.zeros((5, 1)), obs=obs, var=pd.DataFrame(index=["M1"]))
    roi = pd.DataFrame(
        {
            "Id": [1, 2],
            "Name": ["Alpha", "Beta"],
            "type": ["Rectangle", "Polygon"],
            "all_points": [
                "00,00 10,00 10,10 00,10",
                "05,00 15,00 15,10 05,10",
            ],
        }
    )

    result = addROI_omero(
        adata,
        roi,
        subset="one",
        label="roi",
        buffer_roi=4,
        buffer_regions="Alpha",
        n_jobs=1,
        verbose=True,
    )

    assert result.obs.loc["alpha", "roi"] == "Alpha"
    assert result.obs.loc["overlap", "roi"] == "Alpha_Beta"
    assert result.obs.loc["border", "roi"] == "Alpha_border"
    assert result.obs.loc["beta", "roi"] == "Beta"
    assert pd.isna(result.obs.loc["other-image", "roi"])
    assert "Buffering enabled" in capsys.readouterr().out


def test_spatial_interaction_covers_summary_per_image_and_validation(monkeypatch):
    module = importlib.import_module("scimap.plotting.spatial_interaction")

    adata = _plotting_adata()
    adata.uns["spatial_interaction"] = _interaction_data()
    monkeypatch.setattr(module.plt, "show", Mock())
    clustermap = Mock(return_value=Mock())
    monkeypatch.setattr(module.sns, "clustermap", clustermap)

    summary = module.spatial_interaction(
        adata,
        binary_view=True,
        subset_phenotype="A",
        subset_neighbour_phenotype=["B", "A"],
        return_data=True,
    )
    assert summary.index.tolist() == ["A"]
    assert summary.columns.tolist() == ["B", "A"]
    assert np.isnan(summary.loc["A", "B"])
    assert summary.loc["A", "A"] == 1

    per_image = module.spatial_interaction(
        adata,
        summarize_plot=False,
        subset_phenotype=["B", "A"],
        subset_neighbour_phenotype="A",
        p_val=0.5,
        return_data=True,
    )
    assert per_image.index.get_level_values("phenotype").tolist() == ["B", "A"]
    assert per_image.columns.tolist() == ["one", "two"]
    assert clustermap.call_count == 2

    missing = _plotting_adata()
    with pytest.raises(ValueError, match="not found"):
        module.spatial_interaction(missing)
    adata.uns["single"] = _interaction_data().drop(columns=["two", "pvalue_two"])
    with pytest.raises(ValueError, match="single image"):
        module.spatial_interaction(
            adata, spatial_interaction="single", summarize_plot=False
        )


@pytest.mark.parametrize("image_order", [["image_B", "image_A"], ["image_A", "image_B"]])
@pytest.mark.parametrize("pvalue_order", [["image_B", "image_A"], ["image_A", "image_B"]])
def test_spatial_interaction_aligns_pvalues_by_image(
    monkeypatch, image_order, pvalue_order
):
    module = importlib.import_module("scimap.plotting.spatial_interaction")
    data = pd.DataFrame({
        "phenotype": ["T_cell", "B_cell"],
        "neighbour_phenotype": ["Macrophage", "Macrophage"],
        "image_B": [0.8, -0.5],
        "image_A": [1.5, 0.6],
        "pvalue_image_B": [0.9, 0.01],
        "pvalue_image_A": [0.001, 0.9],
    })
    data = data[["phenotype", "neighbour_phenotype"] + image_order
                + [f"pvalue_{name}" for name in pvalue_order]]
    adata = AnnData(np.zeros((2, 1)))
    adata.uns["spatial_interaction"] = data.copy()
    monkeypatch.setattr(module.plt, "show", Mock())
    clustermap = Mock(return_value=Mock())
    monkeypatch.setattr(module.sns, "clustermap", clustermap)

    result = module.spatial_interaction(
        adata, summarize_plot=False, p_val=0.05,
        row_cluster=False, col_cluster=False, return_data=True,
    )

    expected = pd.DataFrame(
        {"image_A": [np.nan, 1.5], "image_B": [-0.5, np.nan]},
        index=pd.MultiIndex.from_tuples(
            [("B_cell", "Macrophage"), ("T_cell", "Macrophage")],
            names=["phenotype", "neighbour_phenotype"],
        ),
    )
    expected = expected.reindex(result.index)
    pd.testing.assert_frame_equal(result, expected)
    pd.testing.assert_frame_equal(clustermap.call_args.kwargs["mask"], expected.isna())
    pd.testing.assert_frame_equal(adata.uns["spatial_interaction"], data)


def test_spatial_distance_covers_heatmap_numeric_distribution_and_validation(
    monkeypatch,
):
    module = importlib.import_module("scimap.plotting.spatial_distance")

    adata = _plotting_adata()
    adata.uns["spatial_distance"] = _distance_data(adata)
    monkeypatch.setattr(module.plt, "show", Mock())
    monkeypatch.setattr(module.sns, "clustermap", Mock(return_value=Mock()))
    catplot = Mock(return_value=Mock())
    displot = Mock(return_value=Mock())
    monkeypatch.setattr(module.sns, "catplot", catplot)
    monkeypatch.setattr(module.sns, "displot", displot)

    heatmap = module.spatial_distance(
        adata,
        log=True,
        subset_col="imageid",
        subset_value="two",
        return_data=True,
    )
    assert heatmap.index.tolist() == ["A", "B", "C"]

    non_summary = module.spatial_distance(
        adata, heatmap_summarize=False, return_data=True
    )
    assert "one_A" in non_summary.index

    numeric = module.spatial_distance(
        adata,
        method="numeric",
        distance_from="A",
        distance_to=["B", "C"],
        x_axis="group",
        y_axis="distance",
        facet_by="imageid",
        plot_type="box",
        return_data=True,
    )
    assert numeric["group"].cat.categories.tolist() == ["B", "C"]
    assert catplot.call_count == 1

    distribution = module.spatial_distance(
        adata,
        method="distribution",
        distance_from="B",
        distance_to="A",
        return_data=True,
    )
    assert distribution["group"].tolist() == ["A", "A"]
    assert displot.call_count == 1

    with pytest.raises(ValueError, match="distance_from"):
        module.spatial_distance(adata, method="numeric", distance_to="A")
    with pytest.raises(ValueError, match="not found"):
        module.spatial_distance(_plotting_adata())
