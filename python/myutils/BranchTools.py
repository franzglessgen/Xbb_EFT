#!/usr/bin/env python
import array

# if a branch does not exist, add it to the tree with a default value
class DefaultIfNotExisting(object):
    
    def __init__(self, branchName, branchType='f', defaultValue=0.0):
        self.debug = False
        self.branches = []
        self.branchName = branchName
        self.branchType = branchType
        self.defaultValue = defaultValue

    def customInit(self, initVars):
        self.sample = initVars['sample']
        self.tree = initVars['tree']
        if not self.tree.GetListOfBranches().FindObject(self.branchName):
            self.branches.append({'name': self.branchName, 'formula': lambda x: self.defaultValue, 'type': self.branchType})

    def getBranches(self):
        return self.branches

class TreeFormula(object):

    def __init__(self, branchName, formula, branchType='f'):
        self.branchName = branchName
        self.branchType = branchType
        self.formula = formula

    def customInit(self, initVars):
        self.sampleTree = initVars['sampleTree']
        self.sampleTree.addFormula(self.formula)
        self.branches = [{'name': self.branchName, 'formula': lambda x: self.sampleTree.evaluate(self.formula), 'type': self.branchType}]

    def getBranches(self):
        return self.branches

class TreeFormulas(object):

    def __init__(self, formulaDict):
        self.formulas = []
        for k,v in formulaDict.iteritems():
            if type(v) == str:
                self.formulas.append(TreeFormula(k, v))
            else:
                self.formulas.append(TreeFormula(k, v['formula'], branchType=v['type']))

    def customInit(self, initVars):
        for x in self.formulas:
            x.customInit(initVars)
    
    def getBranches(self):
        return sum([x.getBranches() for x in self.formulas],[])

class newCutBranches(object):
    def __init__(self, formulaDict):
        self.formulas = []
        self.formulaDict = formulaDict

    def customInit(self, initVars):
        self.config = initVars['config']
        for k,v in self.formulaDict.iteritems():
            print 'k,v', k,v
            if v == 'dc':
                self.formulas.append(TreeFormula(k, self.config.get('Cuts',self.config.get(v+':'+k, 'cut'))))
            #print 'v is', v
            else:
                self.formulas.append(TreeFormula(k, self.config.get(v, k)))
        for x in self.formulas:
            x.customInit(initVars)
    
    def getBranches(self):
        return sum([x.getBranches() for x in self.formulas],[])


class Copy(object):

    def __init__(self, branchName, formula, branchType='f'):
        self.branchName = branchName
        self.branchType = branchType
        self.formula = formula

    def customInit(self, initVars):
        self.branches = [{'name': self.branchName, 'formula': lambda x: getattr(x, self.formula), 'type': self.branchType}]

    def getBranches(self):
        return self.branches

class Collection(object):

    def __init__(self, name, properties, maxSize=1, leaves=False):
        self.name = name
        self.properties = properties
        self.maxSize = maxSize
        self.size = maxSize
        self.scalar = (self.maxSize < 2)
        self.branches = []
        self.branchBuffers = {}
        self.leaves = leaves

        if not self.scalar and not leaves:
            numBranchName = 'n' + self.name
            self.branchBuffers['n'] = array.array('i', [0])
            self.branches.append({'name': numBranchName, 'formula': self.getBranch, 'arguments': 'n', 'type': 'i'})

        # if set to True, the properties are saved as leaves of a branch with name: self.name
        #     warning: at the moment only float supported with leaves
        if leaves:
            # call the central value 'Nominal' in this case
            self.properties = [propertyName if len(propertyName)>0 else 'Nominal' for propertyName in properties]
            self.branchBuffers[name] = array.array('f', [0.0] * len(self.properties))
            leaflist = ':'.join(self.properties) + '/F'
            self.branches.append({
                    'name': name,
                    'formula': self.fillVectorBranch,
                    'arguments': {'branch': name, 'size': len(properties)},
                    'length': len(properties),
                    'leaflist': leaflist,
                    'arrayStyle': True,
                })
            # create index of properties to speed up setting a single property by name
            self.propertyIndex = {}
            for i, propertyName in enumerate(self.properties):
                self.propertyIndex[propertyName] = i

            # add alias for nominal value
            if 'Nominal' in self.propertyIndex and '' not in self.propertyIndex:
                self.propertyIndex[''] = self.propertyIndex['Nominal']

        # otherwise the properties are saved as separate branches with name: self.name + '_' + propertyName
        else:
            for prop in self.properties:
                # properties given as list oif dicts: [{'name':'bla', 'type':'f'},...]
                if type(prop) == dict:
                    propBranch = self.name + '_' + prop['name']
                    propType = prop['type'] if 'type' in prop else 'f'
                    bufferName = prop['name']
                # properties given as list of strings ['bla', ...]
                else:
                    propType = 'f'
                    if prop != '':
                        propBranch = self.name + '_' + prop
                        bufferName = prop
                    else:
                        propBranch = self.name
                        bufferName = ''
                self.branchBuffers[bufferName] = array.array(propType, [0.0] * self.maxSize)
                leaflist = '{branch}[{size}]/{type}'.format(branch=propBranch, size=numBranchName, type=propType.upper()) if not self.scalar else '{branch}/{type}'.format(branch=propBranch, type=propType.upper())
                if self.scalar:
                    self.branches.append({
                            'name': propBranch,
                            'formula': self.fillVectorBranch,
                            'arguments': {'branch': bufferName, 'size': 1},
                            'type': propType,
                            'length': 1,
                            'arrayStyle': True,
                        })
                else:
                    self.branches.append({
                            'name': propBranch,
                            'formula': self.fillVectorBranch,
                            'arguments': {'branch': bufferName, 'size': 'n'},
                            'length': self.maxSize,
                            'leaflist': leaflist,
                        })

    def _b(self, n=None):
        return self.branchBuffers[n] if n else self.branchBuffers[self.name]

    # direct access to branch arrays
    def __getitem__(self, key):
        return self.branchBuffers[key]

    # helper function to set properties which works with leaves=True/False
    #   direct access might be faster for leaves=True
    def setProperty(self, key, value):
        if self.leaves:
            self[self.name][self.propertyIndex[key]] = value
        else:
            self[key][0] = value

    def getProperties(self):
        return self.properties

    def fromList(self, listOfDicts):
        self.setSize(len(listOfDicts))
        for i in range(len(listOfDicts)):
            for p in self.getProperties():
                self[p][i] = listOfDicts[i][p]

    # only for vectors, don't use with scalars
    def setSize(self, size):
        self.branchBuffers['n'][0] = size

    def getSize(self):
        return self.branchBuffers['n'][0]

    def getBranches(self):
        return self.branches

    def getBranch(self, event, arguments):
        return self.branchBuffers[arguments][0]

    def fillVectorBranch(self, event, arguments=None, destinationArray=None):
        size = 1
        if 'size' in arguments:
            if arguments['size'] in self.branchBuffers:
                size = self.branchBuffers[arguments['size']][0]
            elif type(arguments['size']) == int:
                size = arguments['size']
            elif type(arguments['size']) == str:
                size = getattr(event, arguments['size'])
        destinationArray[:size] = self.branchBuffers[arguments['branch']][:size]

class AddCollectionsModule(object):

    def __init__(self):
        self.version = 0
        self.lastEntry = -1
        self.branches = []
        self.branchBuffers = {}
        self.collections = {}
        self.xbbStats = {}

    def customInit(self, initVars):
        for k,v in initVars.items():
            setattr(self, k, v)

    def _r(self, branchName):
        return self.sampleTree._r(branchName)

    def count(self, quantity, increment=1.0):
        if quantity not in self.xbbStats:
            self.xbbStats[quantity] = 0 if type(increment) == int else 0.0
        self.xbbStats[quantity] += increment
   
    # this is called after the last event has been processed, but before the output files are flushed/closed.
    def afterProcessing(self):
        if len(self.xbbStats.keys()) > 0:
            print "INFO: statistics:"
            for k in sorted(self.xbbStats.keys()):
                print " ",k,": ", self.xbbStats[k]

    def getVersion(self):
        return self.version

    def addCollection(self, collection):
        if collection.name not in self.collections:
            self.collections[collection.name] = collection
        else:
            raise Exception("CollectionAlreadyExists")
        self.branches += collection.getBranches()

    def addBranch(self, branchName, default=0.0):
        self.branchBuffers[branchName] = array.array('d', [default])
        self.branches.append({'name': branchName, 'formula': self.getBranch, 'arguments': branchName})
        #print 'self.getBranch', self.getBranch
        #print 'self.branchBuffers[branchName]', self.branchBuffers[branchName]
        #print 'self.branches', self.branches

    def addIntegerBranch(self, branchName, default=0):
        self.branchBuffers[branchName] = array.array('i', [default])
        self.branches.append({'name': branchName, 'formula': self.getBranch, 'arguments': branchName, 'type': 'i'})

    def addVectorBranch(self, branchName, default=0, branchType='d', length=1, leaflist=None):
        self.branchBuffers[branchName] = array.array(branchType, [default]*length)
        newBranch = {
                'name': branchName,
                'formula': self.fillVectorBranch,
                'type': branchType,
                'arguments': {'branch': branchName, 'size': length},
                'length': length,
                'arrayStyle': True,
            }
        if leaflist:
            newBranch['leaflist'] = leaflist
        self.branches.append(newBranch)

    def addIntegerVectorBranch(self, branchName, default=0, length=1):
        self.addVectorBranch(branchName, default=default, length=length, branchType='i')

    def getBranches(self):
        return self.branches

    def getBranch(self, event, arguments):
        self.processEvent(event)
        return self.branchBuffers[arguments][0]

    def setBranch(self, branchName, value):
        self.branchBuffers[branchName][0] = value

    def _b(self, branchName):
        return self.branchBuffers[branchName]
   
    def _isnominal(self, syst):
        return syst is None or syst.lower() == 'nominal'

    def _variations(self, syst=None):
        return ['Up','Down'] if not self._isnominal(syst) else [None]

    def _v(self, n, syst=None, UD=None, separator="_"):
        if self._isnominal(syst): 
            return n
        else:
            return n + "_{syst}".format(syst=syst)+separator+"{UD}".format(UD=UD)

    def fillVectorBranch(self, event, arguments=None, destinationArray=None):
        size = 1
        if 'size' in arguments:
            if type(arguments['size']) == int:
                size = arguments['size']
            elif arguments['size'] in self.branchBuffers:
                size = self.branchBuffers[arguments['size']][0]
            elif type(arguments['size']) == str:
                size = getattr(event, arguments['size'])
        destinationArray[:size] = self.branchBuffers[arguments['branch']][:size]

    def hasBeenProcessed(self, tree):
        return tree.GetReadEntry() == self.lastEntry

    def markProcessed(self, tree):
        self.lastEntry = tree.GetReadEntry()
