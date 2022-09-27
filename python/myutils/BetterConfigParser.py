import os,re,ConfigParser

ConfigParser.SafeConfigParser.readOld = ConfigParser.SafeConfigParser.read

#def ConfigParser.SafeConfigParser.read(optionstr):
#    print optionstr
#    ConfigParser.SafeConfigParser.read(optionstr)

class BetterConfigParser(ConfigParser.SafeConfigParser):

    # Workaround for python readline bug. Details in
    # https://bugzilla.redhat.com/show_bug.cgi?id=1040009
    if not os.isatty(1):
      os.environ['TERM'] = 'dumb'

    # allow string interpolation with environment variables
    def __init__(self, recursiveReplace=True):
        self.recursiveReplace = recursiveReplace
        ConfigParser.SafeConfigParser.__init__(self, os.environ)

    def get(self, section, option, raw=False):
        result = ConfigParser.SafeConfigParser.get(self, section, option, raw=False)
        if not raw:
            result = self.__replaceSectionwideTemplates(result, section=section)
        return result

    def optionxform(self, optionstr):
        '''
        enable case sensitive options in .ini files
        '''
        return optionstr

    def __replaceSectionwideTemplates(self, data, section=None):
        '''
        replace <section|option> with get(section,option) recursivly
        '''
        # allow recursive replacement in the expression itself, e.g. settings = <!Test|settings_<!Test|method!>!>
        if self.recursiveReplace:
            result = None
            findExpression = re.compile("((.*)\<!(.*)\|(.*?)\!>(.*))*")
            while result != data:
                groups = findExpression.search(data).groups()
                if not groups == (None, None, None, None, None): # expression not matched
                    # use . to refer to the current section
                    if groups[2] == "." and section:
                        groups = groups[:2] + (section, ) + groups[3:] 
                    data = self.__replaceSectionwideTemplates(groups[1], section=section) + self.get(groups[2], groups[3]) + self.__replaceSectionwideTemplates(groups[4], section=section)
                else:
                    result = data
        # old behavior
        else:
            result = data
            findExpression = re.compile("((.*)\<!(.*)\|(.*)\!>(.*))*")
            groups = findExpression.search(data).groups()
            if not groups == (None, None, None, None, None): # expression not matched
                result = self.__replaceSectionwideTemplates(groups[1])
                result += self.get(groups[2], groups[3])
                result += self.__replaceSectionwideTemplates(groups[4])

        return result

#    def read(self, optionstr):
#        print 
#        print "Reading configuration in ",optionstr
#        print 
#        config = ConfigParser.SafeConfigParser()
#        return config.read('./${energy}config/paths')

