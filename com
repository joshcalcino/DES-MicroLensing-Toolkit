
#kernprof -l -v -c ‘import parameters; star = parameters.star()’
import parameters; star = parameters.star()
