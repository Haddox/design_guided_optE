def makessm(c="HHH_rd4_0636",value='stabilityscore',center=1,annot=True):
    wtseq = df.query('name == "%s.pdb"' % c)['sequence'].values[0]
    wt=df.query('name == "%s.pdb"' % c)[value].values[0]
    aas='QENHDRKTSAGMLVIWYFP'
    
    out=np.zeros((43,19))
    out[:] = wt
    
    subdf=df.query('cat == "%s"' % c)
    for name, v in zip(subdf['name'], subdf[value]):
        if name.count('_') == 3:
            pos = int(name.split('_')[-1][1:-1])
            aa = name[-1]
            out[pos-1, aas.index(aa)] = v
    
    plt.figure(figsize=(20,7))
    
    if annot:
        annotation = np.where(out > 1, out, '').T
    else:
        annotation=False
    
    if center=='wt': center=wt
    
    hm = sns.heatmap(out.T, xticklabels=['%s%s' % (aa, pos) for aa, pos in zip(wtseq, range(1,44))], yticklabels=[x for x in aas],
               cmap='bwr',center=center,annot=annotation,fmt='.3s')
    plt.plot([0,43],[5,5],color='black')
    plt.plot([0,43],[7,7],color='black')
    plt.plot([0,43],[11,11],color='black')
    plt.plot([0,43],[15,15],color='black')
    plt.plot([0,43],[18,18],color='black')
    
    
    
    for i in range(43):
        plt.scatter(i+0.5,aas.index(wtseq[i])+0.5, s=30,color='black')
        
    cax = plt.gcf().axes[-1]
    cax.plot([0,2],[wt,wt],color='black',linewidth=3)    
    return hm
out = makessm()
plt.title('HHH_rd4_0636.pdb Stability Score',fontsize=15)
plt.tight_layout()
plt.savefig('HHH_rd4_0636.pdb.pdf')