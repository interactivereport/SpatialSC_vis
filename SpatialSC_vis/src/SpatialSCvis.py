import numpy as np
import scanpy as sc
import pandas as pd

import random

def read_and_qc(sample_name, path=None):
    """ This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param path: path to data
    """

    adata = sc.read_visium(path + str(sample_name),
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'

    return adata




def extract_coord_file(adata_segment,ad_sp,sample_id):

    '''This function takes the assignment matrix between cell_id and spot_id from spatial data as input
     and find the assigned spot for each cell. Then assign the spatial coordinates for each cell based on
     the spot coordinates or cell segmented results in each spot.

    :param adata_segment: segmentation results from Tangram function tg.deconvolve_cell_annotations()
    :param ad_sp: spatial data
    :param sample_id: Name of the sample

    '''
    
    adata_segment.obs['spot']=adata_segment.obs['centroids']
    
    for k in range(0,adata_segment.obs.shape[0]):
        adata_segment.obs['spot'][k]=adata_segment.obs['spot'][k][:-2]
        
    ##initialize the coord file
    df = pd.DataFrame(columns = ['cell_id', 'y', 'x'])
    mapping=ad_sp.obsm["tangram_ct_count"]
    
    for i in range(4,mapping.shape[1]):
        
        name=mapping.columns[i]
        if len(np.where(mapping[name]==1)[0])==0: #for cells that are filtered out
            x='Null'
            y='Null'
        
        else:
            spot=mapping.index[np.where(mapping[name]==1)[0][0]]
            spot_mat=adata_segment.obs[adata_segment.obs['spot']==spot]
            
            if spot_mat.shape[0]==0: #for spots only have one cell, so we do not segment on it
                                    #assign the spot coordinates to the cells
                x=float(mapping[mapping.index==spot]['x'][0])
                y=float(mapping[mapping.index==spot]['y'][0])
                
            else: #for spots segmented into multiple cells, randomly assign one of the spatial coordinates for the cell
                n=random.sample(range(0,spot_mat.shape[0]),1)[0]
                x=adata_segment.obs[adata_segment.obs['spot']==spot]['x'][n]
                y=adata_segment.obs[adata_segment.obs['spot']==spot]['y'][n]
        
        df = df.append({'cell_id' : name, 'y' : y, 'x' : x}, ignore_index = True)
    
    df['capture_id']=sample_id
    
    return df



def transfer_coord_sn(name, adata_sn, adata_vis, coord_file_path,flip=False):
    
    
    '''This function is designd to transfer the spatial coordinates csv file from Tangram to the x,y columns in 
    single-nucleus dataset (adata_sn); Also transfer the H&E image previously only in spatial data (adata_vis) to the 
    single-nucleus dataset (adata_sn)

    :param name: name for sample
    :param adata_sn: single-nucleus/single-cell dataset to mapping the spatial coordinates
    :param adata_vis: Visium spatial data
    :param coord_file_path: path saved from Tangram with the spatial coordinates csv file
    :param flip: bool, default is False, set to true, if the H&E image to overlap is flipped on the x-axis


    '''


    ##crreate a spatial key to save spatial image in single-nucleus data

    adata_sn.uns['spatial'][str('spatial_'+name)]=[]
    
    ##initialize the spatial coordinates in single-nucleus data
    adata_sn.obs['y']=np.nan
    adata_sn.obs['x']=np.nan
    adata_sn.obs['capture_id']='nan'
    
    coord_mat=np.zeros((adata_sn.obs.shape[0],2))
    
    df=pd.read_csv(str(coord_file_path)) #read in saved coordinates file from Tangram
    df.index=df['cell_id']
    df=df[df['cell_id'].isin(adata_sn.obs.index)]  #only keep the cell id in the mapping reference file

    
    if df.shape[0]==0:
        raise ValueError("Cell ids are mismatching between two datasets")
    
    mapping_cells=list(df[df['y']!='Null']['cell_id']) #only keep the cells that are not filtered out in Tangram
    
    for f in range(0,len(mapping_cells)):
        adata_sn.obs['y'][str(mapping_cells[f])] = float(df['y'][str(mapping_cells[f])])
        adata_sn.obs['x'][str(mapping_cells[f])] = float(df['x'][str(mapping_cells[f])])
        adata_sn.obs['capture_id'][str(mapping_cells[f])] = df['capture_id'][str(mapping_cells[f])]
 
    if np.sum(adata_sn.obs['y'].isna()==False)!=len(mapping_cells):
        raise ValueError('Not all the mapping cells receive their mapping coordinates')
    
    if np.sum(adata_sn.obs['x'].isna()==False)!=len(mapping_cells):
        raise ValueError('Not all the mapping cells receive their mapping coordinates')
    
    if np.sum(adata_sn.obs['capture_id']!='nan')!=len(mapping_cells):
        raise ValueError('Not all the mapping cells receive their capture id')
    
    coord_mat[:,0]=adata_sn.obs['y']
    if flip==True:
        coord_mat[:,1]=0-adata_sn.obs['x']
    else:
        coord_mat[:,1]=adata_sn.obs['x']
    
    adata_sn.obsm[str('X_spatial_'+name)]=coord_mat
    adata_sn.obsm[str('X_spatial_'+name)][pd.isna(adata_sn.obsm[str('X_spatial_'+name)])]=float(0)

    adata_sn.uns['spatial'][str('spatial_'+name)] = adata_vis.uns['spatial'][str(name)]
    
    sf=adata_sn.uns['spatial'][str('spatial_'+name)]['scalefactors']['tissue_hires_scalef']
    adata_sn.obsm[str('X_spatial_'+name)]=adata_sn.obsm[str('X_spatial_'+name)]/sf
    

    return adata_sn



