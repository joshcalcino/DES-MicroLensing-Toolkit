#kernprof -l -v 'import cutter; file = open("mag_averages.txt", "r"); numbers = cutter.cutter(file)'
#import cutter; file = open("mag_averages.txt", "r"); numbers = cutter.cutter(file)
import data_practice as dp
data = dp.data_practice()
time = data.avg_mag()
