$project_dir = ENV['URA_PROJECT']
$fwk_dir = ENV['URA']+'/AnalysisFW'
$jobid = ENV['jobid']
$batch = ENV.fetch('batch', "0")
$no_batch_submit = ENV.fetch('no_batch_submit', "0")

tools = "#{$fwk_dir}/rake/tools.rb"
require tools

require 'json'

def psub(target, sub)
  ret = proc {|name| 
    puts name
    puts target
    puts sub
    puts name.sub(target, sub)
    name.sub(target, sub)
  }
  return ret
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
  sh "python ctag_scripts/make_ctag_postfit.py #{args.wp}"
end

task :ctag_scan, [:wp] do |t, args|
  Rake::Task["plots/#{$jobid}/ctageff/mass_discriminant/#{args.wp}/MultiDimScan.root"].invoke()
end

$wroking_points = ['csvLoose', 'csvMedium', 'csvTight', 'ctagLoose', 'ctagMedium', 'ctagTight', 'cmvaLoose', 'cmvaMedium', 'cmvaTight']
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
  sh 'python ctag_scripts/CTagEffPlotter.py  --shapes --wps="ctag*" --noPOIpropagation --noLightFit'
  sh 'python ctag_scripts/CTagEffPlotter.py --plots  --shapes --wps="notag" --noPOIpropagation'
  sh 'python ctag_scripts/CTagEffPlotter.py  --shapes --wps="csv*"  --noLightFit --noPOIpropagation '
end

task :ctag_plotfit do |t|
  Rake::Task['ctag_shapes'].invoke()
  Rake::Task['ctag_fitall'].invoke()
  #Rake::Task['breakdown_all'].invoke()
  sh 'python ctag_scripts/make_ctag_tables.py'
  Rake::Task['breakdown_all'].invoke()
  sh "write_csv.py ctag"
  #sh "mv ctag.csv plots/#{jobid}/ctageff/."
end
