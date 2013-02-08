definition = """
// wordlist allowed for parameter's types

default                      #a simple parameter
boolean                      #can be true or false
text                         #free text

list                         #an item list
+ assembly                   #a list of assembly

numeric                      #a numeric value
+ int                        #an integer
+ float                      #a float

file                         #a simple file
+ txt                        #a text file

+ image                      #a simple graphic
+ + png                      #a png
+ + pdf                      #a pdf

+ track                      #a file describing genomic files
+ + bed                      #a bed file
+ + wig                      #a wig file
+ + bam                      #a BAM file
+ + bw                       #a BigWig file
+ + sql                      #a sql file

+ userfile                   #a file that always need to be inputed by th user
assembly                     #a list of assemblies
"""

FILE = 'file'
TRACK = 'track'

wordlist = {}
inclusions = {}
parent_types = []


for line in definition.split('\n'):
    if not line.find('#') == -1:
        word, comment = line.split('#')
        word_stripped, comment = ''.join(c for c in word if c not in '+').strip(), comment.strip()
        wordlist[word_stripped] = comment
        if not word.count('+') > 0:                         # we have a root element
            parents = [word_stripped]
            parent_types.append(word_stripped)
            count = 0

        else:
            c = word.count('+')
            if c > count:                                   # we have a deeper child element

                cur_parent = parents[-1]
                if cur_parent not in inclusions:
                    inclusions[cur_parent] = []
                inclusions[cur_parent].append(word_stripped)
                parents.append(word_stripped)
                count += 1

            elif c == count:                                 # we have a child for same previous parent
                if len(parents) > count:
                    parents.pop()
                cur_parent = parents[-1]
                inclusions[cur_parent].append(word_stripped)
                parents.append(word_stripped)

            else:                                            # (c < count) we go back one step
                parents.pop()
                parents.pop()
                cur_parent = parents[-1]
                if cur_parent not in inclusions:
                    inclusions[cur_parent] = []
                inclusions[cur_parent].append(word_stripped)
                parents.append(word_stripped)
                count -= 1


def is_of_type(obj, oftype):
    """
    Look if obj is of type "oftype":
    """
    if obj == oftype:
        return True
    if oftype not in inclusions:
        return False
    types = inclusions.get(oftype)
    if obj in types:
        return True
    for ty in types:
        if is_of_type(obj, ty):
            return True
    return False


def parent_type(t):
    """
    Get the "parent type of t. The parent type is one of the "top types".
    """
    if t in parent_types:
        return t
    for k, v in inclusions.iteritems():
        if t in v:
            if k in parent_types:
                return k
            return parent_type(k)
    raise Exception("No parent found for type %s" % t)
