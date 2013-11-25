from ivy import treegraph as tg

merged = {}
with open('ncbi/merged.dmp') as f:
    for line in f:
        v = line.split()
        merged[int(v[0])] = int(v[2])

g = tg.load_taxonomy_graph('ncbi/ncbi.xml.gz')
## for v in g.vertices():
##     if g.vertex_name[v]=='Viridiplantae': break
## plants = v
## g = tg.load_taxonomy_graph('ncbi/viridiplantae.filt.xml.gz')

def check_trees():
    outf = open('pb184.filt','w')
    failed = open('pb184.nonparsable','w')
    with open('pb.dmp.nr.trees.184') as f:
        for line in f:
            pbtree, s = line.split()
            try:
                r = tg.ivy.tree.read(s)
                for lf in r.leaves():
                    w = lf.label.split('_')
                    assert w[-1].startswith('ti')
                    lf.taxid = int(w[-1][2:])
                    assert lf.taxid
                for n in r.postiter():
                    n.parent = None; del n.children
                outf.write(line)
            except:
                print 'failed:', pbtree
                failed.write(line)
    outf.close()
    failed.close()

def color_vertices(g, tg, taxv):
    nxt, bck = g.hindex[g.taxid_vertex[taxv]]
    colored = tg.new_vertex_property('bool')
    seen = set()

    lvs = set()
    for v in tg.vertices():
        if v.out_degree() == 1:
            seen.add(v)
            taxv = g.vertex(tg.vertex_taxv[v])
            if g.incertae_sedis[taxv]:
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
    for v in tg.vertices():
        if colored[v]:
            i = check(v)
            if i: c[i].append(v)

    return colored, c

def proc(line, probs):
    pbtree, s = line.split()
    print 'processing', pbtree
    r = tg.ivy.tree.read(s)
    lvs = r.leaves()
    rps = []
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
        probs.write('%s\n' % pbtree)
        return []

    r.mrca = tg.rootpath_mrca(rps)
    r.taxv = g.taxid_vertex[r.mrca]
    taxids = set()
    for rp in rps:
        while 1:
            if rp[-1] == r.mrca: break
            else: rp.pop()
        assert rp
        taxids.update(rp)
    taxidsubg = tg.taxid_subgraph(g, taxids)
    for x in [ taxidsubg.vertex(int(lf.taxv)) for lf in lvs
               if leaf_tid_counts[lf.taxid]==1 ]:
        taxidsubg.vfilt[x] = 0
    
    tg = tg.gt.Graph(directed=False)
    tg.mrca = r.mrca
    print 'mrca:', g.taxid_name(r.mrca)
    tg.vertex_taxid = tg.get_or_create_vp(tg, 'taxid', 'int')
    tg.vertex_taxv = tg.get_or_create_vp(tg, 'taxv', 'int')
    v2lf = {}
    N = len(r)
    verts = tg.add_vertex(N)
    for n in r:
        n.v = verts.next()
        if not n.children:
            tg.vertex_taxid[n.v] = n.taxid
            tg.vertex_taxv[n.v] = int(n.taxv)
            v2lf[n.v] = n
        if n.parent:
            tg.add_edge(n.parent.v, n.v)

    convex = {}

    def traverse(taxv):
        tid = g.vertex_taxid[taxv]
        p, c = color_vertices(g, tg, tid)
        if len(c)==1 and len(c[1])==1: # taxv/tid is convex
            rv = c[1][0] # rv is the root of the convex subtree
            tg.set_vertex_filter(p)
            lvs = [ x for x in tg.vertices() if x.out_degree()==1 ]
            if len(lvs) > 2:
                mrca = tg.rootpath_mrca(
                    [ tg.taxid_rootpath(taxidsubg, tg.vertex_taxid[lf])
                      for lf in lvs ]
                    )
                ancv = [taxidsubg.taxid_vertex[mrca]]
                while ancv[-1] != taxv:
                    ancv.append(ancv[-1].in_neighbours().next())
                k = '.'.join([ str(taxidsubg.vertex_taxid[x])
                               for x in ancv ])
                convex[k] = (rv, p)
            tg.set_vertex_filter(None)
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
        tg.set_vertex_filter(p)
        s = make_newick(root, set([root]))
        tg.set_vertex_filter(None)
        names = ','.join([ g.taxid_name(int(x)) for x in k.split('.') ])
        newicks.append((k, names, s))
        print 'subtree:', names

    for n in r.postiter():
        n.parent = None; del n.children

    return newicks

outf = open('pb184.filt.subtrees','w')
probs = open('pb184.filt.problem_subtrees','w')
with open('pb184.filt') as f:
    ## flag = True
    ## while flag:
    ##     line = f.readline()
    ##     if line.startswith('ti85819_ci81'): flag = False
    for line in f:
        pbtree = line.split()[0]
        for k, names, s in proc(line, probs):
            outf.write('%s\t%s\t%s\t%s;\n' % (pbtree, k, names, s))
outf.close()
probs.close()
