
def transfer_coord_sn(name, adata_sn, adata_vis, coord_file_path):
    
    adata_sn.obs['y']=np.nan
    adata_sn.obs['x']=np.nan
    adata_sn.obs['capture_id']='nan'
    coord_mat=np.zeros((adata_sn.obs.shape[0],2))
    df=pd.read_csv(str(coord_file_path))
    df.index=df['cell_id']
    df=df[df['cell_id'].isin(adata_sn.obs.index)]

    
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
    coord_mat[:,1]=adata_sn.obs['x']
    
    adata_sn.obsm[str('X_spatial_'+name)]=coord_mat
    adata_sn.obsm[str('X_spatial_'+name)][pd.isna(adata_sn.obsm[str('X_spatial_'+name)])]=float(0)

    adata_sn.uns['spatial'][str('spatial_'+name)] = adata_vis.uns['spatial'][str(name)]
    
    sf=adata_sn.uns['spatial'][str('spatial_'+name)]['scalefactors']['tissue_hires_scalef']
    adata_sn.obsm[str('X_spatial_'+name)]=adata_sn.obsm[str('X_spatial_'+name)]/sf
    

    return adata_sn


