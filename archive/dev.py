# %%
import anndata as ad
import scimap as sm

adata = ad.read_h5ad(
    '/Users/aj/Partners HealthCare Dropbox/Ajit Nirmal/nirmal lab/resources/exemplarData/scimapExampleData/scimapExampleData.h5ad'
)
adata
# %%

# %%
image_path = '/Users/aj/Partners HealthCare Dropbox/Ajit Nirmal/nirmal lab/resources/exemplarData/scimapExampleData/registration/exemplar-001.ome.tif'
image_viewer(
    adata=adata,
    image_path=image_path,
    overlay='phenotype',
    embed_notebook=False,
    point_color='white',
)
# %%

# %%
image_viewer(
    adata=adata,
    image_path=image_path,
    # overlay='phenotype',
    embed_notebook=True,
    # point_color='white',
    backend='vitessce',
)
# %%


# %%
from vitessce import (
    VitessceConfig,
    CoordinationLevel as CL,
    get_initial_coordination_scope_prefix,
)
from os.path import join
import ipywidgets

# %%

# %%
vc = VitessceConfig(schema_version="1.0.16", name="BioMedVis Challenge")
dataset = vc.add_dataset(name="Blood Vessel", uid="bv").add_file(
    url="/Users/aj/Partners HealthCare Dropbox/Ajit Nirmal/nirmal lab/resources/exemplarData/scimapExampleData/registration/exemplar-001.ome.tif",
    file_type="image",
)

# spatial = vc.add_view("spatialBeta", dataset=dataset)
# lc = vc.add_view("layerControllerBeta", dataset=dataset)

# vc.layout(spatial | lc)
vw = vc.widget(js_package_version="3.4.5", remount_on_uid_change=False)
vw
# %%


# %%
import vizarr
import zarr

image_path = '/Users/aj/Partners HealthCare Dropbox/Ajit Nirmal/nirmal lab/resources/exemplarData/scimapExampleData/registration/exemplar-001.ome.tif'
image = tiff.TiffFile(image_path, is_ome=False)  # is_ome=False
store = zarr.open(image.aszarr())  # convert image to Zarr array
viewer = vizarr.Viewer()
viewer.add_image(store)
viewer
# %%
