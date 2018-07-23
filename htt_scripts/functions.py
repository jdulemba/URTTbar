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

