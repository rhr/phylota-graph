import ivy
from ivy import treegraph as tg

def color_vertices(taxonomy, treegraph, taxv):
    """
    taxv: NCBI taxon id
    
    Color the vertices of `treegraph` that are members of taxon `taxv`
    """
    nxt, bck = taxonomy.hindex[taxonomy.taxid_vertex[taxv]]
    colored = treegraph.new_vertex_property('bool')
    seen = set()

    lvs = set()
    for v in treegraph.vertices():
        if v.out_degree() == 1:
            seen.add(v)
            taxv = g.vertex(treegraph.vertex_taxv[v])
            if taxonomy.incertae_sedis[taxv]:
                p = taxv.in_neighbours().next()
                pn, pb = g.hindex[p]
                if nxt >= pn and bck <= pb:
                    colored[v] = 1
                    lvs.add(v)
            else:
                n, b = g.hindex[taxv]
                if n >= nxt and b <= bck:
                    colored[v] = 1
                    lvs.add(v)

    def gather():
        s = set()
        for v in lvs:
            for n in v.out_neighbours():
                if not colored[n]: s.add(n)
        return s

    def check(x):
        i = 0
        for x in x.out_neighbours():
            if not colored[x]: i += 1
        return i

    verts = gather()
    while 1:
        for x in verts:
            if check(x) == 1:
                lvs.add(x)
                colored[x] = 1
        v = gather()
        if v == verts: break
        verts = v

    c = tg.defaultdict(list)
    for v in treegraph.vertices():
        if colored[v]:
            i = check(v)
            if i: c[i].append(v)

    return colored, c

def proc(g, line, merged, probfile, outfile):
    pbtree, s = line.split()
    print 'processing', pbtree
    r = ivy.newick.parse(s)
    lvs = r.leaves()
    rps = [] # rootpaths: list of taxids from leaf to root
    leaf_tid_counts = tg.Counter()
    try:
        for lf in lvs:
            w = lf.label.split('_')
            lf.gi = int(w[-2][2:])
            lf.taxid = int(w[-1][2:])
            leaf_tid_counts[lf.taxid] += 1
            if lf.taxid not in g.taxid_vertex and lf.taxid in merged:
                lf.taxid = merged[lf.taxid]
            lf.taxv = g.taxid_vertex[lf.taxid]
            lf.taxid_next, lf.taxid_back = g.hindex[lf.taxv]
            lf.taxid_rootpath = tg.taxid_rootpath(g, lf.taxid)
            for i, x in enumerate(lf.taxid_rootpath):
                if x not in g.taxid_vertex and x in merged:
                    lf.taxid_rootpath[i] = merged[x]
            rps.append(lf.taxid_rootpath)
    except:
        print '!!! problem assigning leaf taxids'
        probfile.write('%s\n' % pbtree)
        #return []

    r.mrca = tg.rootpath_mrca(rps)
    r.taxv = g.taxid_vertex[r.mrca]
    taxids = set()
    for rp in rps:
        # trim rootpaths: make them terminate with mrca
        while 1:
            if rp[-1] == r.mrca: break
            else: rp.pop()
        assert rp
        taxids.update(rp)
    taxidsubg = tg.taxid_subgraph(g, taxids)

    # no need to check for convexity for singleton tip taxa
    for x in [ taxidsubg.vertex(int(lf.taxv)) for lf in lvs
               if leaf_tid_counts[lf.taxid]==1 ]:
        taxidsubg.vfilt[x] = 0
    
    treegraph = tg.gt.Graph(directed=False)
    treegraph.mrca = r.mrca
    print 'mrca:', g.taxid_name(r.mrca)
    treegraph.vertex_taxid = tg.get_or_create_vp(treegraph, 'taxid', 'int')
    treegraph.vertex_taxv = tg.get_or_create_vp(treegraph, 'taxv', 'int')
    v2lf = {}
    N = len(r)
    verts = treegraph.add_vertex(N)
    for n in r:
        n.v = verts.next()
        if not n.children:
            treegraph.vertex_taxid[n.v] = n.taxid
            treegraph.vertex_taxv[n.v] = int(n.taxv)
            v2lf[n.v] = n
        if n.parent:
            treegraph.add_edge(n.parent.v, n.v)

    convex = {}

    def traverse(taxv):
        tid = g.vertex_taxid[taxv]
        p, c = color_vertices(g, treegraph, tid)
        if len(c)==1 and len(c[1])==1: # taxv/tid is convex
            rv = c[1][0] # rv is the root of the convex subtree
            treegraph.set_vertex_filter(p)
            lvs = [ x for x in treegraph.vertices() if x.out_degree()==1 ]
            if len(lvs) > 2:
                mrca = tg.rootpath_mrca(
                    [ tg.taxid_rootpath(taxidsubg, treegraph.vertex_taxid[lf])
                      for lf in lvs ]
                    )
                ancv = [taxidsubg.taxid_vertex[mrca]]
                while ancv[-1] != taxv:
                    # STRANGE EDGE CASES HERE
                    try: ancv.append(ancv[-1].in_neighbours().next())
                    except StopIteration: pass

                k = '.'.join([ str(taxidsubg.vertex_taxid[x])
                               for x in ancv ])
                convex[k] = (rv, p)
            treegraph.set_vertex_filter(None)
        else:
            for n in taxv.out_neighbours():
                traverse(n)

    for v in taxidsubg.root.out_neighbours(): traverse(v)

    def make_newick(root, seen):
        children = [ x for x in root.out_neighbours() if x not in seen ]
        if children:
            seen.update(children)
            s = '(%s)' % ','.join(
                [ make_newick(c, seen) for c in children ]
                )
        else:
            s = v2lf[root].label.replace(',','').replace('(','').replace(')','')
        return s
        
    newicks = []

    for k, (root, p) in convex.items():
        treegraph.set_vertex_filter(p)
        s = make_newick(root, set([root]))
        treegraph.set_vertex_filter(None)
        names = ','.join([ g.taxid_name(int(x)) for x in k.split('.') ])
        outfile.write('%s\t%s\t%s\t%s;\n' % (pbtree, k, names, s))
        print 'wrote subtree:', names

    for n in r.postiter():
        n.parent = None; del n.children


if __name__ == "__main__":
    merged = {}
    with open('ncbi/merged.dmp') as f:
        for line in f:
            v = line.split()
            merged[int(v[0])] = int(v[2])

    g = tg.load_taxonomy_graph('ncbi/ncbi.xml.gz')

    probfile = open('pb/pb184.readable.problem_subtrees','w')
    outfile = open('pb/pb184.readable.convex_subtrees','w')
    with open('pb/pb184.readable.trees') as f:
        for line in f:
            proc(g, line, merged, probfile, outfile)
    outfile.close()
    probfile.close()
