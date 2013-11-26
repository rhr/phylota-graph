# Some trees from phylota have tip labels that corrupt the newick
# string. This script filters them out.

import ivy, re

outf = open('pb/pb184.readable.trees','w')
failed = open('pb/pb184.unreadable.trees','w')
with open('pb/pb.dmp.maximalnr.trees.184') as f:
    for line in f:
        for bad, good in (('_(','_'),(')_','_'),(',_','_'),('BOLD:','BOLD_')):
            line = line.replace(bad, good)
        pbtree, s = line.split()
        try:
            r = ivy.newick.parse(s)
            for lf in r.leaves():
                w = lf.label.split('_')
                assert w[-1].startswith('ti')
                lf.taxid = int(w[-1][2:])
                assert lf.taxid
            assert r.parent is None
            for n in r.postiter():
                n.parent = None; del n.children
            outf.write(line)
        except:
            print 'failed:', pbtree
            failed.write(line)
outf.close()
failed.close()
