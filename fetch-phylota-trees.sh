#!/bin/bash

src=http://phylota.net/pb/Download
data=pb.dmp.maximalnr.trees.184.gz
dest=pb
mkdir $dest
cd $dest
wget $src/$data
gunzip $data
