def extract_coord_file(adata_segment,ad_sp):
    
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

