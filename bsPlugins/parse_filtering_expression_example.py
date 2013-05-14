
import re

def parse_filter(filter):
    def operation(f):
        try:
            if not f: return lambda x:True
            fsplit = re.split('(or|and|and not|or not)',f)
            if re.search('[><!=]',fsplit[0]):
                expression = ''
                for j,e in enumerate(fsplit):
                    if j%2==0: # operation
                        expression += 'float(x)'+e
                    else: # boolean operator
                        expression += ' %s '%e
                return lambda x: eval(expression)
            else:
                return lambda x: bool(re.search(f,x))
        except:
            raise ValueError("Unknown expression: %s" % f)
    filter = filter.split(',')
    for k,f in enumerate(filter):
        filter[k] = operation(f)
    return filter

def parse_splitline(splitline,filter):
    """Return a tuple of booleans of the same length as the input tuple,
    according to whether each value passed its filter or not."""
    return tuple(filter[k](x) for k,x in enumerate(splitline))

###################################  TEST  ###################################

def test():
    test = "../tests/testing_files/genes_table.tab"
    filterline = "ENSMUS,>301 and < 1000, >200 and <10000"
    filter = parse_filter(filterline)
    f = open(test,'r')
    f.readline()
    for line in f:
        splitline = line.split()
        print parse_splitline(splitline,filter)
    f.close()

test()

