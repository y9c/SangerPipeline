#!/usr/bin/env python
# encoding: utf-8
import re

map_table = {}
for l in open("./hehe_table", "r").readlines()[1:]:
    ll = re.split('\W+', l)
    map_table[ll[0]] = ll[1]
    print ll

print map_table
