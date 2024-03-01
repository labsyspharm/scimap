#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Fri Feb 26 19:47:10 2021
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.tl.spatial_lda`: This function constructs a neighborhood matrix based on 
    user-specified categorical variables, such as cell types, 
    and applies Latent Dirichlet Allocation (LDA) to model the latent space of 
    cellular distributions. It returns weights that describe the spatial 
    organization of cells, facilitating the identification of Recurrent Cellular Neighborhoods (RCNs). 
    
    The `sm.tl.spatial_cluster` function should be utilized to cluster these 
    latent vectors into RCNs, offering insights into the spatial dynamics 
    of cellular environments.

## Function
"""

#Import
from sklearn.neighbors import BallTree
import numpy as np
import pandas as pd
import re

# Gensim
import gensim
import gensim.corpora as corpora
from gensim.models import CoherenceModel

# Function
def spatial_lda (adata, 
                 x_coordinate='X_centroid',
                 y_coordinate='Y_centroid',
                 z_coordinate= None,
                 phenotype='phenotype', 
                 method='radius', 
                 radius=30, 
                 knn=10,
                 imageid='imageid',
                 num_motifs=10, 
                 random_state=0, 
                 subset=None,
                 verbose=True,
                 label='spatial_lda',**kwargs):
    """
Parameters:
        adata (anndata.AnnData):  
            AnnData object, containing spatial gene expression data.

        x_coordinate (str, required):  
            Column name in `adata` denoting the x-coordinates.

        y_coordinate (str, required):  
            Column name in `adata` denoting the y-coordinates.

        z_coordinate (str, optional):  
            Column name in `adata` for z-coordinates, for 3D spatial data.

        phenotype (str, required):  
            Column name in `adata` indicating cell phenotype or classification.

        method (str, optional):  
            Neighborhood definition method: 'radius' for fixed distance, 'knn' for K nearest neighbors.

        radius (int, optional):  
            Radius defining local neighborhoods (when method='radius').

        knn (int, optional):  
            Number of nearest neighbors for neighborhood definition (when method='knn').

        imageid (str, optional):  
            Column name in `adata` specifying image identifiers, for analyses within specific images.

        num_motifs (int, optional):  
            Number of latent motifs to identify.

        random_state (int, optional):  
            Seed for random number generator, ensuring reproducibility.

        subset (str, optional):  
            Specific image identifier for targeted analysis.
        
        verbose (bool, optional):  
            If True, enables progress and informational messages.

        label (str, optional):  
            Custom label for storing results in `adata.uns`.

Returns:
        adata (anndata.AnnData):  
            The input `adata` object, updated with spatial LDA results in `adata.uns[label]`.

Example:
        ```python
        
        # Analyze spatial motifs using the radius method 
        adata = sm.tl.spatial_lda(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                            method='radius', radius=50, num_motifs=10,
                            label='lda_radius_50')
    
        # KNN method with specific image subset
        adata = sm.tl.spatial_lda(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                            method='knn', knn=15, num_motifs=15, subset='image_01',
                            label='lda_knn_15_image_01')
    
        # 3D spatial data analysis using the radius method
        adata = am.tl.spatial_lda(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid', z_coordinate='Z_centroid',
                            method='radius', radius=100, num_motifs=20, label='lda_3D_radius_100')
        ```
    
    """
    

    # Function
    def spatial_lda_internal (adata_subset, x_coordinate,y_coordinate,z_coordinate,phenotype, 
                              method, radius, knn, imageid):
        
        # Print which image is being processed
        if verbose:
            print('Processing: ' + str(np.unique(adata_subset.obs[imageid])))
        
        
        # Create a dataFrame with the necessary inforamtion
        if z_coordinate is not None:
            if verbose:
                print("Including Z -axis")
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'z': adata_subset.obs[z_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        else:
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})

        
        # Create a DataFrame with the necessary inforamtion
        #data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        
        # Identify neighbourhoods based on the method used
        # a) KNN method
        
        if method == 'knn':
            if verbose:
                print("Identifying the " + str(knn) + " nearest neighbours for every cell")
            if z_coordinate is not None:
                tree = BallTree(data[['x','y','z']], leaf_size= 2)
                ind = tree.query(data[['x','y','z']], k=knn, return_distance= False)
            else:
                tree = BallTree(data[['x','y']], leaf_size= 2)
                ind = tree.query(data[['x','y']], k=knn, return_distance= False)
            ind = list(np.array(item) for item in ind)
                

        # b) Local radius method
        if method == 'radius':
            if verbose:
                print("Identifying neighbours within " + str(radius) + " pixels of every cell")
            if z_coordinate is not None:
                kdt = BallTree(data[['x','y','z']], metric='euclidean') 
                ind = kdt.query_radius(data[['x','y','z']], r=radius, return_distance=False)
            else:
                kdt = BallTree(data[['x','y']], metric='euclidean') 
                ind = kdt.query_radius(data[['x','y']], r=radius, return_distance=False)


# =============================================================================
#         if method == 'knn':
#             if verbose:
#                 print("Identifying the " + str(knn) + " nearest neighbours for every cell")
#             tree = BallTree(data[['x','y']], leaf_size= 2)
#             ind = tree.query(data[['x','y']], k=knn, return_distance= False)
#             #ind = [np.array(x) for x in ind]
#             ind = list(np.array(item) for item in ind)
#             
#         # b) Local radius method
#         if method == 'radius':
#             if verbose:
#                 print("Identifying neighbours within " + str(radius) + " pixels of every cell")
#             kdt = BallTree(data[['x','y']], leaf_size= 2) 
#             ind = kdt.query_radius(data[['x','y']], r=radius, return_distance=False)
#             
# =============================================================================

        # Map phenotype
        phenomap = dict(zip(list(range(len(ind))), data['phenotype'])) # Used for mapping
        for i in range(len(ind)):
            ind[i] = [phenomap[letter] for letter in ind[i]]
            
        # return
        return ind
    
    # Subset a particular image if needed
    if subset is not None:
        adata_list = [adata[adata.obs[imageid] == subset]]
    else:
        adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
        
    # Apply function to all images
    # Create lamda function 
    r_spatial_lda_internal = lambda x: spatial_lda_internal(adata_subset=x,
                                                            x_coordinate=x_coordinate,
                                                            y_coordinate=y_coordinate,
                                                            z_coordinate=z_coordinate,
                                                            phenotype=phenotype, 
                                                            method=method, 
                                                            radius=radius, 
                                                            knn=knn, 
                                                            imageid=imageid) 
    all_data = list(map(r_spatial_lda_internal, adata_list)) # Apply function 
    
    # combine all the data into one
    texts = np.concatenate( all_data, axis=0 ).tolist()
    
    # LDA pre-processing
    if verbose:
        print ('Pre-Processing Spatial LDA')
    # Create Dictionary
    id2word = corpora.Dictionary(texts)

    # Term Document Frequency
    corpus = [id2word.doc2bow(text) for text in texts]
    
    # Build LDA model
    if verbose:
        print ('Training Spatial LDA')
    try:
        lda_model = gensim.models.ldamulticore.LdaMulticore(corpus=corpus,
                                                   id2word=id2word,
                                                   num_topics=num_motifs, 
                                                   random_state=random_state,**kwargs)
    except:
        lda_model = gensim.models.ldamodel.LdaModel(corpus=corpus,
                                                   id2word=id2word,
                                                   num_topics=num_motifs, 
                                                   random_state=random_state,**kwargs)
    
    # Compute Coherence Score
    if verbose:
        print ('Calculating the Coherence Score')
    coherence_model_lda = CoherenceModel(model=lda_model, texts=texts, dictionary=id2word, coherence='c_v')
    coherence_lda = coherence_model_lda.get_coherence()
    if verbose:
        print('\nCoherence Score: ', coherence_lda)

    # isolate the latent features
    if verbose:
        print ('Gathering the latent weights')
    topic_weights = []
    for row_list in lda_model[corpus]:
        tmp = np.zeros(num_motifs)
        for i, w in row_list:
            tmp[i] = w
        topic_weights.append(tmp)
    # conver to dataframe
    arr = pd.DataFrame(topic_weights, index=adata.obs.index).fillna(0)
    arr = arr.add_prefix('Motif_')
    
    # isolate the weights of phenotypes
    pattern = "(\d\.\d+).\"(.*?)\""
    cell_weight = pd.DataFrame(index=np.unique(adata.obs[phenotype]))
    for i in range(0, len(lda_model.print_topics())):
        level1 = lda_model.print_topics()[i][1]
        tmp = pd.DataFrame(re.findall(pattern, level1))
        tmp.index = tmp[1]
        tmp = tmp.drop(columns=1)
        tmp.columns = ['Motif_'+ str(i)]
        cell_weight = cell_weight.merge(tmp, how='outer', left_index=True, right_index=True)
    # fill zeros
    cell_weight = cell_weight.fillna(0).astype(float)
    
    # save the results in anndata object
    adata.uns[label] = arr # save the weight for each cell
    adata.uns[str(label)+'_probability'] = cell_weight # weights of each cell type
    #adata.uns[str(label)+'_model'] = lda_model
    
    # return
    return adata
