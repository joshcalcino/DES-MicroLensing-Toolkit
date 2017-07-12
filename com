<<<<<<< HEAD
#kernprof -l -v 'import cutter; file = open("mag_averages.txt", "r"); numbers = cutter.cutter(file)'
#import cutter; file = open("mag_averages.txt", "r"); numbers = cutter.cutter(file)
import data_practice as dp
data = dp.data_practice()
time = data.avg_mag()
=======

#kernprof -l -v -c ‘import parameters; star = parameters.star()’
import parameters; star = parameters.star()
>>>>>>> 8accf0a2de7404b28571f65e61c47bf7d5e9fd6b
