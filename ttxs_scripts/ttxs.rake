$project_dir = ENV['URA_PROJECT']
$fwk_dir = ENV['URA']+'/AnalysisFW'
$jobid = ENV['jobid']
$batch = ENV.fetch('batch', "0")
$no_batch_submit = ENV.fetch('no_batch_submit', "0")

tools = "#{$fwk_dir}/rake/tools.rb"
require tools

require 'json'

def psub(target, sub)
  return proc {|name| name.sub(target, sub)}
end

def psubs(target, *subs)
  return subs.map{|x| psub(target, x)}
end

task :new_ttbar_trial, [:info] do |t, args|
  new_trial('ttxs', 'ttxsec', args.info)
end

task :new_ttbar_plots, [:info] do |t, args|
  new_trial('', 'ttxsec', args.info)
end

task :publish_ttxsec do |t|
  link = `ls -ld plots/#{$jobid}/ttxsec`.scan(/-> (.+)/).flatten.last
  if not link
    link = 'ttxsec'
  end
  publish_pics("plots/#{$jobid}/#{link}", "#{ENV['HOME']}/public_html/#{$jobid}/#{link}")
end

rule /\.model\.root$/ => psub(/\.model\.root$/, '.txt') do |t|
  dir = File.dirname(t.name)
  json_file = File.basename(t.source).sub(/.txt$/, '.json')
  chdir(dir) do
    sh "text2workspace.py #{File.basename(t.source)} --X-allow-no-background -P URAnalysis.AnalysisTools.statistics.TTBarXsecFitter:ttxsecfitterWJetCategories --PO yieldsJson=#{json_file} -o #{File.basename(t.name)}"
  end
end

$external_toys=''
rule /(:?\.toy)?\.mlfit\.root$/ => psub(/(:?\.toy)?\.mlfit\.root$/, '.model.root') do |t|
  dir = File.dirname(t.name)
  toy_cmd = ''
  if t.name.include? '.toy.'
    toy_cmd = '--saveToys --expectSignal 1 -t 200'
  else
    #if there are no toys make plots!
    toy_cmd = '--saveShapes'
  end

  if not $external_toys.empty?
    toy_cmd += " --toysFile #{$external_toys}"
  end

  chdir(dir) do
    combine_cmd = "combine #{File.basename(t.source)} -M MaxLikelihoodFit --saveNormalizations --saveWithUncertainties --saveNLL --skipBOnlyFit --minos=all"
    if $batch == "1" and t.name.include? '.toy.'
      sh "rm -f mlfit[0-9]*.root higgsCombine[0-9]*.MaxLikelihoodFit.mH120.[0-9]*.root"
      sh "cp #{ENV['URA']}/AnalysisTools/scripts/batch_job.sh ."
      File.open('condor.mltoys.jdl', 'w') do |file|
        file << "universe = vanilla\n"
        file << "Executable = batch_job.sh\n"
        file << "Should_Transfer_Files = YES\n"
        file << "WhenToTransferOutput = ON_EXIT\n"
        file << "Transfer_Input_Files = #{File.basename(t.source)}\n"
        
        (2345678...2745678).step(10000) do |seed|
          file << "\n"
          file << "Output = con_#{seed}.stdout\n"
          file << "Error = con_#{seed}.stderr\n"
          file << "Arguments = #{combine_cmd} #{toy_cmd} -n #{seed} -s #{seed}\n"
          file << "Queue\n"
        end
      end
      if $no_batch_submit == "0"
        sh 'condor_submit condor.mltoys.jdl'
        sh 'hold.py --check_correctness=./ --maxResubmission=0'
        sh "merge_toys.py #{File.basename(t.name)} mlfit[0-9]*.root" 
        sh "merge_toy_experiments.py toy_experiments.root higgsCombine[0-9]*.MaxLikelihoodFit.mH120.[0-9]*.root"
      end
    else
      sh "#{combine_cmd} #{toy_cmd} &> #{File.basename(t.name).sub('.root','.log')}"
      if t.name.include? '.toy.'
        sh "mv higgsCombineTest.MaxLikelihoodFit.mH120.123456.root toy_experiments.root"
      end
      sh "cp mlfit.root #{File.basename(t.name)}"
    end
  end
end

rule /\.harvested\.root$/ => psub(/\.harvested\.root$/, '.mlfit.root') do |t|
  toy_cmd = ''
  if t.name.include? '.toy.'
    toy_cmd = '-t'
  end
  sh "./harvest_fit.py #{t.source} #{t.source.sub(/(:?\.toy)?\.mlfit\.root$/, '.binning.json')} #{t.source.sub(/(:?\.toy)?\.mlfit\.root$/, '.root')} -o #{t.name} #{toy_cmd}"
end

def get_toy_harvest(fname)
  harvesdir = fname.split('toys')[0]
  var = fname.sub("plots/#{$jobid}/ttxsec/", "").split('/')[0]
  return "#{harvesdir}/#{var}.toy.harvested.root"
end

$quick_toy_check='' #urgh, global vars, really? Yes...
$bias_cmd=''
rule /toys\/shapes\/tt_right.json$/ => proc {|name| get_toy_harvest(name)} do |t|
  vardir = File.dirname(t.source)
  var = File.basename(t.source).sub('.toy.harvested.root','')
  sh "./toy_diagnostics.py #{var} #{t.source} #{vardir}/#{var}.toy.mlfit.root -o #{vardir}/toys #{$quick_toy_check} #{$bias_cmd}"
end

task :toy_diagnostics, [:vardir, :biasID] do |t, args|
  puts args.vardir, args.biasID
  bias_potfix=''
  if args.biasID 
    bias_harvest="plots/#{$jobid}/ttxsec_#{args.biasID}/#{args.vardir}/#{args.vardir}.toy.harvested.root"
    Rake::Task[bias_harvest].invoke
    $bias_cmd="--biasFile=#{ENV['URA_PROJECT']}/#{bias_harvest}"
    bias_potfix+="_#{args.biasID}"
  end
  
  Rake::Task["plots/#{$jobid}/ttxsec/#{args.vardir}#{bias_potfix}/toys/shapes/tt_right.json"].invoke
end

task :fit, [:var] do |t, args|
  Rake::Task["plots/#{$jobid}/ttxsec/#{args.var}/#{args.var}.harvested.root"].invoke
end

task :fit_all do |t|
  candidates = Dir.glob("plots/#{$jobid}/ttxsec/*/*.txt")
  vars = []
  candidates.each do |candidate|
    txt_var = File.basename(candidate, '.txt')
    dir_var = File.basename(File.dirname(candidate))
    if txt_var == dir_var
      vars << txt_var
    end
  end
  
  vars.each do |var|
    Rake::Task["plots/#{$jobid}/ttxsec/#{var}/#{var}.harvested.root"].invoke
  end
end

task :toys, [:var] do |t, args|
  Rake::Task["plots/#{$jobid}/ttxsec/#{args.var}/#{args.var}.toy.mlfit.root"].invoke
end

task :fit_toys, [:var] do |t, args|
  Rake::Task["plots/#{$jobid}/ttxsec/#{args.var}/#{args.var}.toy.harvested.root"].invoke
end

task :fit_toys_wbias, [:var, :biasID] do |t, args|
  #Make toys with bias (provided in a different trial directory)
  Rake::Task["plots/#{$jobid}/ttxsec_#{args.biasID}/#{args.var}/#{args.var}.toy.mlfit.root"].invoke
  
  #create a clone of the variable directory
  #with the datacard and model (if any)
  new_varname = "#{args.var}_#{args.biasID}"
  sh "mkdir -p plots/#{$jobid}/ttxsec/#{new_varname}"
  chdir("plots/#{$jobid}/ttxsec/#{new_varname}") do
    sh "cp ../#{args.var}/#{args.var}.txt #{new_varname}.txt"
    sh "cp ../#{args.var}/#{args.var}.root #{new_varname}.root"
    sh "cp ../#{args.var}/#{args.var}.binning.json #{new_varname}.binning.json"
    sh "sed -i 's|#{args.var}|#{new_varname}|g' #{new_varname}.binning.json"
    if File.exists? "../#{args.var}/#{args.var}.model.root"
      sh "cp ../#{args.var}/#{args.var}.model.root #{new_varname}.model.root"
    end
  end
  
  sh "python make_skewed_matrices.py #{args.var} #{args.biasID}"

  #disable batch, we are quick
  #set external toy source
  biased_toys = "plots/#{$jobid}/ttxsec_#{args.biasID}/#{args.var}/toy_experiments.root"
  puts biased_toys
  if not File.exists? biased_toys
    throw "The biased toys do not exist!"
  end
  $external_toys="#{ENV['URA_PROJECT']}/#{biased_toys}"
  pbatch = $batch
  $batch="0"

  #run fit
  Rake::Task["plots/#{$jobid}/ttxsec/#{new_varname}/#{new_varname}.toy.harvested.root"].invoke

  #clean up the environment
  $external_toys=''
  $batch=pbatch
end

task :postfit_plots, [:var] do |t, args|
  sh "python make_postfit_plots.py #{args.var}"
end

task :postfit_all do |t|
  candidates = Dir.glob("plots/#{$jobid}/ttxsec/*/*.txt")
  vars = []
  candidates.each do |candidate|
    txt_var = File.basename(candidate, '.txt')
    dir_var = File.basename(File.dirname(candidate))
    if txt_var == dir_var
      vars << txt_var
    end
  end
  
  vars.each do |var|
    Rake::Task["plots/#{$jobid}/ttxsec/#{var}/#{var}.harvested.root"].invoke
    sh "python make_postfit_plots.py #{var}"
  end
end

task :prefit_plots, [:var] do |t, args|
  sh "python make_prefit_plots.py #{args.var}"
end

task :prefit_all do |t|
  candidates = Dir.glob("plots/#{$jobid}/ttxsec/*/*.txt")
  vars = []
  candidates.each do |candidate|
    txt_var = File.basename(candidate, '.txt')
    dir_var = File.basename(File.dirname(candidate))
    if txt_var == dir_var
      vars << txt_var
    end
  end
  
  vars.each do |var|
    sh "python make_prefit_plots.py #{var}"
  end
end


#
#   BIN OPTIMIZATION
#

def optimize_bin(var, varmin, vmin, vmax, vstep, vrange, bin_name='Bin0', predecessors='')
  $quick_toy_check = '--nopars --no-post-pulls'
  prev_batch = $batch
  $batch = '1'
  prev_no_batch_submit = $no_batch_submit
  $no_batch_submit = "1"
  varmin = Float(varmin)
  max_edges = []
  (Float(vmin)..Float(vmax)).step(Float(vstep)) do |varmax|
    max_edges << varmax
  end
  sh "python runTTXSecPlotter.py --optimize_binning '#{var}:[#{predecessors}]:#{varmin}:[#{max_edges.join(',')}]' --subdir=binning_optimization"

  bin_dirs = []
  (Float(vmin)..Float(vmax)).step(Float(vstep)) do |varmax|
    vardir = "#{predecessors.gsub(',','_')}#{format("%.1f", varmin)}_#{format("%.1f", varmax)}"
    bin_dir = "plots/#{$jobid}/ttxsec/#{var}/binning_optimization/#{vardir}"
    mlfits = "#{bin_dir}/#{var}.toy.mlfit.root"
    #creates model, creates condor.jdl DOES NOT SUBMIT IT
    Rake::Task[mlfits].invoke
    bin_dirs << bin_dir
  end

  bin_dirs.each do |bin_dir|
    chdir(bin_dir) do
      sh 'condor_submit condor.mltoys.jdl' 
    end
  end

  bin_dirs.each do |bin_dir|
    chdir(bin_dir) do  
      sh "hold.py --check_correctness=./ --maxResubmission=0"
      sh "merge_toys.py #{var}.toy.mlfit.root mlfit[0-9]*.root"
    end
  end

  jsons = []
  bin_dirs.each do |bin_dir|
    jfile = "#{bin_dir}/toys/shapes/tt_right.json"
    sh "date"
    Rake::Task[jfile].invoke
    jsons << jfile
  end
  opt_dir = "plots/#{$jobid}/ttxsec/#{var}/binning_optimization"
  sh "./plot_optimization.py #{var} #{opt_dir} #{vrange} #{jsons.join(' ')} --name=#{bin_name}"
  $quick_toy_check=''
  $no_batch_submit=prev_no_batch_submit
  $batch = prev_batch
  return "#{opt_dir}/#{bin_name}.json"
end

task :optimize_binning, [:var, :varmin, :vmin, :vmax, :vstep] do |t, args|
  optimize_bin(args.var, args.varmin, args.vmin, args.vmax, args.vstep, Float(args.vmax)-Float(args.vmin))
end

task :optimize_bins, [:var, :vmin, :vmax, :vstep, :maxsize, :already_done] do |t, args|
  low_bounds = []
  if args.already_done
    low_bounds = args.already_done.split('_').map{|x| Float(x)}
  end
  low_bounds << Float(args.vmin)
  vmax = Float(args.vmax)
  maxsize = Float(args.maxsize)
  vstep = Float(args.vstep)
  var = args.var
  while low_bounds[-1] < vmax do
    limit = low_bounds[-1]+maxsize < vmax ? low_bounds[-1]+maxsize : vmax 
    preced = ""
    if low_bounds.size > 1
      preced = "#{low_bounds[0..-2].join(',')},"
    end
    sh "date"
    sh "echo starting optimization of bin #{low_bounds.size}" 
    json_file = optimize_bin(
                             var, 
                             low_bounds[-1], 
                             low_bounds[-1]+vstep, 
                             limit, 
                             vstep, 
                             maxsize,
                             "Bin#{low_bounds.size}", 
                             preced
      )
    jmap = JSON.parse(File.open(json_file).read)
    low_bounds << jmap['best']
  end
  sh "echo #{low_bounds.join(' , ')} > plots/#{$jobid}/ttxsec/#{var}/#{var}.optimalbin.txt"
end

=begin
def optimise_nbins(steering_file)
  #TODO: write steering
  
  niter = 0
  maxiter = 500
  while niter < maxiter do
    niter += 1
    jmap = JSON.parse(File.open(steering_file).read)
    binning = jmap['binning']
    var = jmap['var']
    opt_dir = "plots/#{$jobid}/ttxsec/#{var}/binning_optimization"
    sh "echo #{steering_file} >> #{opt_dir}/iterations.txt"
    output_dir = "#{opt_dir}/#{binning.map{|x| format("%.1f", x)}.join('_')}"
    output_txt = "#{output_dir}/#{var}.txt"
    if File.exist?(output_txt)
      puts 'Circular reference found! Exiting'
      break
    end

    sh "python runTTXSecPlotter.py --binning '#{var}:[#{binning.join(',')}]' --subdir=binning_optimization"
    jfile = "#{output_dir}/toys/shapes/tt_right.json"
    Rake::Task[jfile].invoke

    sh "./find_next_iteration.py #{steering_file} #{jfile} #{output_dir}/steering.json"
    steering_file="#{output_dir}/steering.json"
  end
end
=end
