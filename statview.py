#viewing cProfile output from
#python3 -m cProfile -o prof.stats icing.py
import pstats

p = pstats.Stats('prof.stat')

p.sort_stats('time').print_stats(35)
#p.sort_stats('call').print_stats(15)
