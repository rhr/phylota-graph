experimental code for extracting 'convex subtrees' (rooted subtrees
corresponding to a monophyletic taxon) from unrooted phylota trees

requires latest git version of ivy

# fetch the NCBI taxonomy
sh fetch-taxonomy.sh

# make a taxonomy graph
python make-taxonomy-graph.py

# fetch phylota trees (pb184)
sh fetch-phylota-trees.sh

# filter out problematic trees (can't be read due to newick-unfriendly labels)
python filter-trees.py

# extract convex subtrees from readable phylota trees
python extract-convex-subtrees.py
