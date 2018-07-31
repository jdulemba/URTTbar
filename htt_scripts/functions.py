from pdb import set_trace
import rootpy.plotting as plotting

def stack_plots(lists):

    lists.sort(key=lambda x: x.Integral())
    total = 0
    ratio_hists = []
    Stack_hists = []

    for i in lists:
        total += i
    for i in lists:
        ratio_hists.append(i/total)
        #i.SetFillStyle(1001)
        Stack_hists.append(i)

    stack = plotting.HistStack()
    norm_stack = plotting.HistStack()

#    print 'stack plots'
#    set_trace()

    for i in Stack_hists:
        stack.Add(i)
        norm_stack.Add(i/total)


    return stack, norm_stack, ratio_hists

def print_table(lines, filename, separate_head=True):
    """Prints a formatted table given a 2 dimensional array"""
    #Count the column width
    widths = []
    for line in lines:
        for i,size in enumerate([len(x) for x in line]):
            while i >= len(widths):
                widths.append(0)
            if size > widths[i]:
                widths[i] = size

    #Generate the format string to pad the columns
    print_string = ""
    for i,width in enumerate(widths):
        print_string += "{" + str(i) + ":" + str(width) + "} | "
    if (len(print_string) == 0):
        return
    print_string = print_string[:-3]

    with open(filename, 'w') as f:
        #Print the actual data
        for i,line in enumerate(lines):
            print >> f, print_string.format(*line)
            if (i == 0 and separate_head):
                print >> f, "-"*(sum(widths)+3*(len(widths)-1))



