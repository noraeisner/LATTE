import lightkurve as lk

search_result = lk.search_tesscut(searchtic, sector=sec)
tpf = search_result.download(cutout_size=15)

print (search_result)
print (tpf)