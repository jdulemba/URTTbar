$project_dir = ENV['URA_PROJECT']
$fwk_dir = ENV['URA']+'/AnalysisFW'
$jobid = ENV['jobid']
$batch = ENV.fetch('batch', "0")
$no_batch_submit = ENV.fetch('no_batch_submit', "0")

tools = "#{$fwk_dir}/rake/tools.rb"
require tools

analysis_rules = "#{$fwk_dir}/rake/analysis.rake"
import analysis_rules

meta_rules = "#{$fwk_dir}/rake/meta.rake"
import meta_rules

htt_rules = "#{$project_dir}/htt_scripts/htt.rake"
import htt_rules

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

task :set_trial,[:trialID] do |t, args|
  chdir("results/#{$jobid}") do
    link  = "ctag_eff"
    results = "#{link}_#{args.trialID}"
    if File.exist? link
      if not File.symlink? link
        throw "#{link} MUST be a symlink in the first place to work!"
      end
      sh "unlink #{link}"
    end
  
    if not File.exist? results
      throw "trial ID not found!"
    end

    sh "ln -s #{results} #{link}"
  end

  chdir("plots/#{$jobid}") do
    link  = "ctageff"
    results = "#{link}_#{args.trialID}"
    if File.exist? link
      if not File.symlink? link
        throw "#{link} MUST be a symlink in the first place to work!"
      end
      sh "unlink #{link}"
    end
    
    if not File.exist? results
      throw "trial ID not found!"
    end

    sh "ln -s #{results} #{link}"
  end
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

###########################
# Charm-tagging tasks
###########################

task :new_ctag_trial, [:info] do |t, args|
  new_trial('ctag_eff', 'ctageff', args.info)
end

task :new_ctag_plots, [:info] do |t, args|
  new_trial('', 'ctageff', args.info)
end

task :publish_ctag do |t|
  link = `ls -ld plots/#{$jobid}/ctageff`.scan(/-> (.+)/).flatten.last
  if not link
    link = 'ctageff'
  end
  publish_pics("plots/#{$jobid}/#{link}", "#{ENV['HOME']}/public_html/#{$jobid}/#{link}")
end

rule /fitModel.root$/ => psub(/fitModel.root$/, 'datacard.txt') do |t|
  dir = File.dirname(t.name)
  chdir(dir) do
    puts 'creating workspace'
    if File.readlines("datacard.txt").grep(/NOLIGHTSFFIT/).size > 0
      opts = '--PO fitLightEff=False --PO lightConstantsJson=datacard.json'
    else
      opts = '--PO lightConstantsJson=datacard.json'
    end
    if File.readlines("datacard.txt").grep(/NOPOIPROPAGATION/).size > 0
      opts += ' --PO POIPropagation=False'
    end
    sh "text2workspace.py datacard.txt -P URAnalysis.AnalysisTools.statistics.CTagEfficiencies:ctagEfficiency #{opts} -o #{File.basename(t.name)}"
  end
end

task :ctag_model, [:wp] do |t, args|
  Rake::Task["plots/#{$jobid}/ctageff/mass_discriminant/#{args.wp}/fitModel.root"].invoke()
end

rule /MultiDimFit(:?Toy|Asimov)?.root$/ => psub(/MultiDimFit(:?Toy|Asimov)?.root$/, 'fitModel.root') do |t|
  toy_cmd = ''
  seed = ""
  if t.name.include? 'Toy'
    toy_cmd = '--saveToys --expectSignal 1 -t 200'
  elsif t.name.include? 'Asimov'
    toy_cmd = '--saveToys --expectSignal 1 -t -1'
    seed = ".123456"
  end
  dir = File.dirname(t.name)
  chdir(dir) do
    puts 'running Multi-dimensional with Profile-Likelyhood errors'#--setPhysicsModelParameterRanges charmSF=0,2:lightSF=0,2
    sh "combine fitModel.root -M MultiDimFit --algo=singles --saveWorkspace #{toy_cmd} > #{File.basename(t.name).sub(/\.root$/,'.log')}"
    sh "cat #{File.basename(t.name).sub(/\.root$/,'.log')}"
    sh "mv higgsCombineTest.MultiDimFit.mH120#{seed}.root #{File.basename(t.name)}"
  end
end

rule /MultiDimScan(:?Toy|Asimov)?.root$/ => psub(/MultiDimScan(:?Toy|Asimov)?.root$/, 'fitModel.root') do |t|
  toy_cmd = ''
  seed = ''
  if t.name.include? 'Toy'
    toy_cmd = '--saveToys --expectSignal 1 -t 200'
  elsif t.name.include? 'Asimov'
    toy_cmd = '--saveToys --expectSignal 1 -t -1'
    seed = ".123456"
  end
  dir = File.dirname(t.name)
  chdir(dir) do
    puts 'running Multi-dimensional with Profile-Likelyhood errors'
    sh "combine -M MultiDimFit fitModel.root --algo=grid --points=400 --setPhysicsModelParameterRanges charmSF=0.05,2.05:lightSF=0.05,2.05 #{toy_cmd} > #{File.basename(t.name).sub(/\.root$/,'.log')}"
    sh "mv higgsCombineTest.MultiDimFit.mH120#{seed}.root #{File.basename(t.name)}"
  end
end

rule /MaxLikeFit(:?Toy|Asimov)?.root$/ => psub(/MaxLikeFit(:?Toy|Asimov)?.root$/, 'fitModel.root') do |t|
  dir = File.dirname(t.name)
  bname = File.basename(t.name)
  combine_cmd = "combine fitModel.root -M MaxLikelihoodFit --saveNormalizations --saveWithUncertainties --saveWorkspace --minos=all --saveNLL  --skipBOnlyFit"
  chdir(dir) do
    if t.name.include? 'Toy'
      if $batch == "1"
        toy_cmd = '--saveToys --expectSignal 1 -t 10'
        sh "cp #{ENV['URA']}/AnalysisTools/scripts/batch_job.sh ."
        File.open('condor.mltoys.jdl', 'w') do |file|
          file << "universe = vanilla\n"
          file << "Executable = batch_job.sh\n"
          file << "Should_Transfer_Files = YES\n"
          file << "WhenToTransferOutput = ON_EXIT\n"
          file << "Transfer_Input_Files = #{Dir.pwd}/fitModel.root\n"
          
          (2345678...2545678).step(10000) do |seed|
            file << "\n"
            file << "Output = con_#{seed}.stdout\n"
            file << "Error = con_#{seed}.stderr\n"
            file << "Arguments = #{combine_cmd} #{toy_cmd} -n #{seed} -s #{seed}\n"
            file << "Queue\n"
          end
        end
        sh 'condor_submit condor.mltoys.jdl'
        sh 'hold.py --check_correctness=./ --maxResubmission=0'
        sh "merge_toys.py #{bname} mlfit[0-9]*.root" 
      else
        toy_cmd = '--saveToys --expectSignal 1 -t 200 -v -1'
        sh "#{combine_cmd} #{toy_cmd}"
        sh "mv mlfit.root #{File.basename(t.name)}"      
      end
    else
      toy_cmd = '--saveShapes '
      if t.name.include? 'Asimov'
        toy_cmd += '--saveToys --expectSignal 1 -t -1'
      end
      puts 'running MaxLikelihood fit with Profile-Likelyhood errors'
      sh "#{combine_cmd} #{toy_cmd} &> fit.log"
      sh "cat fit.log"
      sh "mv mlfit.root #{File.basename(t.name)}"      
      sh "mv higgsCombineTest.MaxLikelihoodFit.mH120.root MLFit_workspace.root"
    end
    #sh "mv higgsCombineTest.MultiDimFit.mH120.root MultiDimFit.root"
  end
end


rule /MaxLikeFitStatOnly\.root$/ => psub(/MaxLikeFitStatOnly/, 'MultiDimFit') do |t|
  dir = File.dirname(t.name)
  bname = File.basename(t.name)
  chdir(dir) do
    sh "combine MultiDimFit.root -M MaxLikelihoodFit -S 0 --minos=all --snapshotName MultiDimFit"
    sh "mv mlfit.root #{File.basename(t.name)}"      
  end
end

task :sys_breakdown, [:wp] do |t, args|
  wpdir = "plots/#{$jobid}/ctageff/mass_discriminant/#{args.wp}"
  Rake::Task["#{wpdir}/MultiDimFit.root"].invoke()
  Rake::Task["#{wpdir}/MaxLikeFitStatOnly.root"].invoke()
  sh "mkdir -p #{wpdir}/sys_breakdown/"  
  nuisances=File.readlines("#{wpdir}/datacard.txt").grep(/( lnN )|( shape )|( param )/).map {|x| x.split()[0]}
  #groups = [/_bin_/, /_MCStat/]
  #group_names = ['bbb', 'MCStat']
  singles = [/JES/,/pu/]
  chdir(wpdir) do
    sh "combine MultiDimFit.root -M MaxLikelihoodFit -S 0 --minos=all --snapshotName MultiDimFit &> /dev/null"
    sh "mv mlfit.root MaxLikeFitStatistic.root"
    nuisances.each do |nuisance|
      if singles.map {|g| g =~ nuisance}.any?
        puts nuisance
        to_freeze = nuisances.select{|x| x != nuisance}.join(',')
        sh "combine MultiDimFit.root -M MaxLikelihoodFit --minos=all --snapshotName MultiDimFit --freezeNuisances=#{to_freeze} &> /dev/null "
        sh "mv mlfit.root sys_breakdown/#{nuisance}.root"
      end
    end
    
    puts 'running for groups'
    ## groups.zip(group_names).each do |group, name|
    ##   to_freeze = nuisances.select{|x| not (group =~ x)}.join(',')
    ##   if to_freeze.length == 0
    ##     next
    ##   end
    ##   sh "combine MultiDimFit.root -M MaxLikelihoodFit --minos=all --snapshotName MultiDimFit --freezeNuisances=#{to_freeze} &> /dev/null"
    ##   sh "mv mlfit.root sys_breakdown/#{name}.root"
    ## end
    to_freeze = nuisances.select{|x| singles.map{|y| y =~ x}.any? }.join(',')
    sh "combine MultiDimFit.root -M MaxLikelihoodFit --minos=all --snapshotName MultiDimFit --freezeNuisances=#{to_freeze} &> /dev/null"
    sh "mv mlfit.root sys_breakdown/other.root"
  end
end

#tasks
task :ctag_fit, [:wp] do |t, args|
  Rake::Task["plots/#{$jobid}/ctageff/mass_discriminant/#{args.wp}/MaxLikeFit.root"].invoke()
end

task :ctag_postfit, [:wp] do |t, args|
  Rake::Task["plots/#{$jobid}/ctageff/mass_discriminant/#{args.wp}/MaxLikeFit.root"].invoke()
  sh "python make_ctag_postfit.py #{args.wp}"
end

task :ctag_scan, [:wp] do |t, args|
  Rake::Task["plots/#{$jobid}/ctageff/mass_discriminant/#{args.wp}/MultiDimScan.root"].invoke()
end

$wroking_points = ['csvLoose', 'csvMedium', 'csvTight', 'ctagLoose', 'ctagMedium', 'ctagTight']
task :ctag_fitall do |t|
  $wroking_points.each do |wp|
    Rake::Task["ctag_postfit"].invoke(wp)
    Rake::Task["ctag_postfit"].reenable
  end
end

task :breakdown_all do |t|
  $wroking_points.each do |wp|
    Rake::Task['sys_breakdown'].invoke(wp)
    Rake::Task["sys_breakdown"].reenable
  end
end

task :ctag_toys, [:wp] do |t, args|
  Rake::Task["plots/#{$jobid}/ctageff/mass_discriminant/#{args.wp}/MaxLikeFitToy.root"].invoke()
end

task :ctag_toy_diagnostics, [:wp ] do |t, args|
  toy_dir = "plots/#{$jobid}/ctageff/mass_discriminant/#{args.wp}"
  sh "mkdir -p #{toy_dir}/toys"
  sh "python toy_diagnostics.py '' '' #{toy_dir}/MaxLikeFitToy.root -o #{toy_dir}/toys/ --use-prefit --noshapes --filter-out-pars='.*_bin_\d+$'"
end

task :ctag_shapes do |t|
  sh 'python CTagEffPlotter.py  --shapes --wps="ctag*" --noPOIpropagation'
  sh 'python CTagEffPlotter.py --plots  --shapes --wps="notag" --noPOIpropagation'
  sh 'python CTagEffPlotter.py  --shapes --wps="csv*"  --noLightFit --noPOIpropagation '
end

task :ctag_plotfit do |t|
  Rake::Task['ctag_shapes'].invoke()
  Rake::Task['ctag_fitall'].invoke()
  #Rake::Task['breakdown_all'].invoke()
  sh 'python make_ctag_tables.py'
  #sh "write_csv.py ctag"
  #sh "mv ctag.csv plots/#{jobid}/ctageff/."
end
