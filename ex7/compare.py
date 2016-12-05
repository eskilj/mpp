import filecmp

correct = filecmp.cmp('output/out_edge192x128_4.pgm', 'output/ser_out_edge192x128_4.pgm')
print correct
