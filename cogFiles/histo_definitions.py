"""Generates code for EmJetHistoMaker
Conventions:

user_XXX functions are used to define specfic instances of elements.
generate_XXX functions are called from cog to generate code fragments.
All other functions are used internally."""

from histo_functions import *
import array

standalone = False
if __name__=='__main__':
    print "Running standalone mode - Printing to screen"
    standalone = True
else:
    try: import cog
    except: ImportError

def outputline(line):
    if not standalone:
        cog.outl(line)
    else:
        print(line)

############################################################
# Object definitions
############################################################
# namedtuple objects can be used as lightweight classes - similar to passive structs in C++

def user_define_histos():
    """Define histograms in this function"""
    histo_dict = OrderedDict()
    name = 'ht'                        ; histo_dict[name] = Histo1F(name , Bins(100 , 0   , 2500 ) )
    name = 'jet_pt'                    ; histo_dict[name] = Histo1F(name , Bins(150 , 0   , 1500 ) )
    name = 'jet_eta'                   ; histo_dict[name] = Histo1F(name , Bins(100 , -5  ,    5 ) )
    name = 'jet_phi'                   ; histo_dict[name] = Histo1F(name , Bins(100 , -5  ,    5 ) )
    name = 'jet_nTrack'                ; histo_dict[name] = Histo1F(name , Bins(100 ,  0. ,  100 ) )
    name = 'nJet_tag'                  ; histo_dict[name] = Histo1F(name , Bins( 10 ,  0  ,   10 ) )
    name = 'n2tag_0'                   ; histo_dict[name] = Histo1F(name , Bins(300 ,  0  ,   30 ) )
    name = 'n2tag_1'                   ; histo_dict[name] = Histo1F(name , Bins(300 ,  0  ,   30 ) )
    name = 'n2tag_2'                   ; histo_dict[name] = Histo1F(name , Bins(300 ,  0  ,   30 ) )
    name = 'n2tag_3'                   ; histo_dict[name] = Histo1F(name , Bins(300 ,  0  ,   30 ) )
    name = 'n2tag_4'                   ; histo_dict[name] = Histo1F(name , Bins(300 ,  0  ,   30 ) )
    name = 'n2tag_5'                   ; histo_dict[name] = Histo1F(name , Bins(300 ,  0  ,   30 ) )
    name = 'n2tag_6'                   ; histo_dict[name] = Histo1F(name , Bins(300 ,  0  ,   30 ) )

    # Jet plot variations
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if name[:4]=='jet_':
            histo_clone = clone_object(histo, postfix='B')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='L')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)

    for name, histo in histo_dict.iteritems():
        if name[:4]=='jet_':
            histo_clone = clone_object(histo, postfix='Emerging')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='Standard')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)

    for name, histo in histo_dict.iteritems():
        if name[:4]=='jet_' or name[:2]=='ht':
            histo_clone = clone_object(histo, postfix='0tag')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='1tag')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='2tag')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='Predicted0To2Tag')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='QCDPredicted1To2Tag')
            histo_clone_dict[histo_clone.name] = histo_clone
            histo_clone = clone_object(histo, postfix='GJetPredicted1To2Tag')
            histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)

    '''
    histo_clone_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
         #histo_clone = clone_object(histo, postfix='Jetfiltered')
         #histo_clone_dict[histo_clone.name] = histo_clone
    histo_dict.update(histo_clone_dict)    
    '''
    return histo_dict

def user_define_histo_vectors():
    """Define histogram vectors in this function"""
    histo_dict = user_define_histos()
    histo_vector_dict = OrderedDict()
    return histo_vector_dict

def construct_bin_str(bins):
    # Fixed width binning
    if type(bins).__name__ == 'Bins':
        return ', {nBins}, {min}, {max}'.format(nBins=bins.nBins, min=bins.min, max=bins.max)
    # Variable width binning
    if type(bins).__name__ == 'VBins':
        return ', nBins_{binname}, bins_{binname}'.format(nBins=bins.nBins, binname=bins.binname)

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

def generate_histo_map_init():
    for name, histo in user_define_histos().iteritems():
        outputline( construct_histo_map_init(histo) )

def calculate_index(histo_dict):
    """Takes histo_dict and calculates index in vectors based on their order in dictionary and their type (1D vs 2D)"""
    index_1d = 0
    index_2d = 0
    histo_id_dict = OrderedDict()
    for name, histo in histo_dict.iteritems():
        if type(histo).__name__ == 'Histo1F':
            histo_id_dict[name] = index_1d
            index_1d += 1
        elif type(histo).__name__ == 'Histo2F':
            histo_id_dict[name] = index_2d
            index_2d += 1
    return histo_id_dict

def generate_histo_index():
    histo_id_dict = calculate_index( user_define_histos() )
    for name, histo in user_define_histos().iteritems():
        template_str = 'if (name=="{name}") return {index};'
        index = histo_id_dict[name]
        output_str = template_str.format(name=name, index=index)
        outputline(output_str)

if __name__=='__main__':
    print 'generate_histo_index():'
    generate_histo_index()
    print''
    print 'generate_histo_map_init():'
    generate_histo_map_init()
    print''

