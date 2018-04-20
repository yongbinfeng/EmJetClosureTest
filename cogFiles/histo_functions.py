"""Generates code for EmJetHistoMaker
Conventions:

user_XXX functions are used to define specfic instances of elements.
generate_XXX functions are called from cog to generate code fragments.
All other functions are used internally."""

from collections import namedtuple
Bins    = namedtuple('Bins'    , ['nBins' , 'min'      , 'max'     , ]) # Fixed width binning
Histo1F = namedtuple('Histo1F' , ['name'  , 'binsX'     , ])           # binsX or binsY contains objects of type Bins or VBins
Histo2F = namedtuple('Histo2F' , ['name'  , 'binsX'    , 'binsY'   , ])
Histo1D = namedtuple('Histo1D' , ['name'  , 'binsX'     , ])           # binsX or binsY contains objects of type Bins or VBins

from collections import OrderedDict

def pad_str(input_string, length=20):
    '''Pads input string with spaces (on the right hand side) to given length'''
    format_string = '{: <%d}' % length
    output_string = format_string.format(input_string)
    return output_string

from collections import namedtuple

from collections import OrderedDict

def compute_fixed_bins(nBins, xMin, xMax):
    delta = xMax - xMin
    binWidth = delta/nBins
    binLowEdges = [(xMin + binWidth * i) for i in xrange(nBins+1)]
    return binLowEdges

def get_type_str(obj):
    if   type(obj).__name__ == 'Histo1F': return 'TH1F'
    elif type(obj).__name__ == 'Histo2F': return 'TH2F'
    elif type(obj).__name__ == 'Histo1D': return 'TH1D'

def get_map_str(obj):
    if   type(obj).__name__ == 'Histo1F': return 'hist1d'
    elif type(obj).__name__ == 'Histo2F': return 'hist2d'
    if   type(obj).__name__ == 'Histo1D': return 'hist1d_double'

def histo_combine1Dto2D(histoX, histoY):
    """Takes two TH1F objects as input and returns a TH2F by combining the inputs"""
    name = '%s_VS_%s' % (histoY.name, histoX.name)
    newhisto = Histo2F(name,
            binsX = histoX.binsX,
            binsY = histoY.binsX, )
    return newhisto

def clone_object(object, prefix='', postfix='', separator='__'):
    """Clone input namedtuple object that has attribute \"name\", by adding prefix and postfix to object.name.
    Separator is the character used to separate the original name from the prefix/suffix."""
    if (not prefix) and (not postfix):
        print 'clone_object: Warning! Either prefix or postfix should be provided.'
        print object
        return object
    else:
        prefix_str  = ''
        postfix_str  = ''
        if prefix  : prefix_str  = prefix + separator
        if postfix : postfix_str = separator + postfix
    newname = '%s%s%s' % (prefix_str, object.name, postfix_str)
    output = object._replace(name=newname)
    return output

def clone_list(input_list, prefix='', postfix=''):
    """Clone list of objects that have attribute \"name\", by adding prefix and postfix to object.name"""
    output_list = []
    for obj in input_list:
        output = clone_object(obj, prefix, postfix)
        output_list.append(output)
    return output_list

def vectorize_histo(histo, prefixlist=[''], postfixlist=['']):
    """Takes an input histogram object and produces a list of histograms with all given combinations of prefixes/postfixes.
    Prefix/postfix list can be a list of strings or numbers. (Numbers automatically get converted to strings.)"""
    output = []
    for prefix in prefixlist:
        if not isinstance(prefix, str): prefix = str(prefix) # Convert to string
        for postfix in postfixlist:
            if not isinstance(postfix, str): postfix = str(postfix) # Convert to string
            histo_clone = clone_object(histo, prefix=prefix, postfix=postfix)
            output.append(histo_clone)
    return output

def offset_bins(bins):
    """Offsets Bins() object with integer bin edges, so that bin centers fall on integers instead"""
    if bins.nBins != bins.max - bins.min: print "Warning: Trying to offset non-integer bins!"
    return Bins(bins.nBins, bins.min-0.5, bins.max-0.5)

def offset(histo):
    '''Use to offset bins by -0.5 to center bins on integer values'''
    return histo._replace( binsX=offset_bins(histo.binsX) )

def construct_bin_str(bins):
    # Fixed width binning
    if type(bins).__name__ == 'Bins':
        return ', {nBins}, {min}, {max}'.format(nBins=bins.nBins, min=bins.min, max=bins.max)

def construct_histo_init(obj):
    name    = obj.name
    typestr = get_type_str(obj)
    binstr  = construct_bin_str(obj.binsX)
    if hasattr(obj, 'binsY'): binstr += construct_bin_str(obj.binsY)
    if hasattr(obj, 'binsZ'): binstr += construct_bin_str(obj.binsZ)
    output = '{name} = new {typestr}("{name}", "{name}" {binstr});'.format(name=name, typestr=typestr, binstr=binstr)
    return output

def construct_histo_map_init(obj):
    name    = obj.name
    typestr = get_type_str(obj)
    mapstr  = get_map_str(obj)
    binstr  = construct_bin_str(obj.binsX)
    if hasattr(obj, 'binsY'): binstr += construct_bin_str(obj.binsY)
    if hasattr(obj, 'binsZ'): binstr += construct_bin_str(obj.binsZ)
    output = '{mapstr}["{name}"] = new {typestr}("{name}", "{name}" {binstr});'.format(name=name, typestr=typestr, mapstr=mapstr, binstr=binstr)
    return output

def generate_histo_decl():
    """Outputs lines like:
    TH1F* jet_pt;
    """
    for name, histo in user_define_histos().iteritems():
        name    = histo.name
        typestr = get_type_str(histo)
        outputline( '{typestr}* {name};'.format(name=name, typestr=typestr) )

def generate_histo_init():
    for name, histo in user_define_histos().iteritems():
        outputline( construct_histo_init(histo) )

def generate_histo_map_init():
    for name, histo in user_define_histos().iteritems():
        outputline( construct_histo_map_init(histo) )

def generate_histo_dest():
    """Outputs lines like:
    delete jet_pt;
    """
    for name, histo in user_define_histos().iteritems():
        name    = histo.name
        typestr = get_type_str(histo)
        outputline( 'delete {name};'.format(name=name, typestr=typestr) )

def generate_histo_vector_decl():
    """Outputs lines like:
    vector<TH1F*> jet_pt_sorted_by_pt;
    """
    for name, histo_vector in user_define_histo_vectors().iteritems():
        typestr = get_type_str(histo_vector[0])
        outputline( 'vector<{typestr}*> {name};'.format(name=name, typestr=typestr) )

def generate_histo_vector_init():
    for name, histo_vector in user_define_histo_vectors().iteritems():
        outputline('{')
        vector_name    = name
        for histo in histo_vector:
            histo_name = histo.name
            outputline( "  auto " + construct_histo_init(histo) )
            outputline( '  {vector_name}.push_back({histo_name});'.format(histo_name=histo_name, vector_name=vector_name) )
        outputline('}')

def generate_histo_vector_dest():
    for name, histo_vector in user_define_histo_vectors().iteritems():
        outputline( 'for (auto i: {name}) {{ delete i; }}'.format(name=name) )
        outputline( '{name}.clear();'.format(name=name) )

if __name__=='__main__':
    print 'generate_histo_decl():'
    generate_histo_decl()
    print''
    print 'generate_histo_init():'
    generate_histo_init()
    print''
    print 'generate_histo_dest():'
    generate_histo_dest()
    print''
    print 'generate_histo_vector_decl():'
    generate_histo_vector_decl()
    print''
    print 'generate_histo_vector_init():'
    generate_histo_vector_init()
    print''
    print 'generate_histo_vector_dest():'
    generate_histo_vector_dest()
    print''

