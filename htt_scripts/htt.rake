$jobid = ENV['jobid']
$project_dir = ENV['URA_PROJECT']

# This is a test task to get a feel for ruby.
task :output do |t|
    sh "ls #{$project_dir}/bin"
end

task :btag_effs do |t|
    py_path = "#{$project_dir}/htt_scripts"
    sh "python #{py_path}/make_btag_effs.py baselinJ20"
#    command = `python make_btag_effs.py baselineJ20`
#    sh "#{command}"
end

$leptons = ['muons', 'electrons']
$njet_opts = ['3', '4+']
task :htt_plots, [:lep, :njets] do |t, args|
    lep_validity = $leptons.include?(args.lep)
    if lep_validity == false
      puts "#{args.lep} isn't a valid lepton choice. Your options are"
      $leptons.each do |lep|
        puts "   #{lep}"
      end
    end
    njet_validity = $njet_opts.include?(args.njets)
    if njet_validity == false
      puts "#{args.njets} isn't a valid jet choice. Your options are"
      $njet_opts.each do |njets|
        puts "   #{njets}"
      end
    end

    sh "python htt_scripts/HTTPlotter.py #{args.lep} --plots --njets=#{args.njets}"

end

task :htt_presel, [:lep, :njets] do |t, args|
    lep_validity = $leptons.include?(args.lep)
    if lep_validity == false
      puts "#{args.lep} isn't a valid lepton choice. Your options are"
      $leptons.each do |lep|
        puts "   #{lep}"
      end
    end
    njet_validity = $njet_opts.include?(args.njets)
    if njet_validity == false
      puts "#{args.njets} isn't a valid jet choice. Your options are"
      $njet_opts.each do |njets|
        puts "   #{njets}"
      end
    end

    sh "python htt_scripts/HTTPlotter.py #{args.lep} --preselection --njets=#{args.njets}"

end
task :all_htt_plots_presel do |t|
    $leptons.each do |lep|
        $njet_opts.each do |njets|
            Rake::Task["htt_presel"].invoke(lep, njets)
            Rake::Task["htt_presel"].reenable
            Rake::Task["htt_plots"].invoke(lep, njets)
            Rake::Task["htt_plots"].reenable
        end
    end
end


task :publish_htt do |t|
  link = `ls -ld plots/#{$jobid}/htt`.scan(/-> (.+)/).flatten.last
  if not link
    link = 'htt'
  end
  publish_pics("plots/#{$jobid}/#{link}", "#{ENV['HOME']}/public_html/#{$jobid}/#{link}")
end

task :publish_tsgen do |t|
  link = `ls -ld plots/#{$jobid}/tsgen`.scan(/-> (.+)/).flatten.last
  if not link
    link = 'tsgen'
  end
  publish_pics("plots/#{$jobid}/#{link}", "#{ENV['HOME']}/public_html/#{$jobid}/#{link}")
end

task :htt_trial,[:trialID] do |t, args|
  chdir("results/#{$jobid}") do
    link  = "htt_simple"
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
    link  = "htt"
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

def cmb2lep(cmb, replace)
  base = File.dirname(File.dirname(File.expand_path(cmb)))
  return "#{base}/#{replace}"
end

rule /htt\/cmb\/cmb\.root/ => [proc{|x| cmb2lep(x, 'muons/muons.root')}, proc{|x| cmb2lep(x, 'electrons/electrons.root')}] do |t|
  sh "mkdir -p plots/#{$jobid}/htt/cmb"
  sh "hadd -f -O #{t.name} #{t.prerequisites.join(" ")}"
end

def timestamp2input(stamp)
  channel_dir = File.dirname(File.dirname(File.expand_path(stamp)))
  channel = File.basename(channel_dir)
  return "#{channel_dir}/#{channel}.root"
end

$widths = ['5', '10', '25', '50']
rule "httlimit.timestamp" => [proc{|x| timestamp2input(x)}, "#{$project_dir}/htt_scripts/setup_limit_htt.py"] do |t|
  channel_dir = File.dirname(File.dirname(File.expand_path(t.name)))
  chan = File.basename(channel_dir)
  chdir("plots/#{$jobid}/htt/#{chan}") do
    $widths.each do |width|
      sh "#{$project_dir}/htt_scripts/setup_limit_htt.py #{chan}.root #{chan} #{width}"
    end
    chdir("output") do
      $widths.each do |width|
        sh "combineTool.py -M T2W -i pseudoscalar_#{width}pc/* -o workspace.root -P CombineHarvester.CombineTools.InterferenceModel:interferenceModel"
      end
      sh 'combineTool.py -M Asymptotic -d */*/workspace.root --there -n .limit --parallel 8 --run blind'
      sh 'combineTool.py -M CollectLimits */*/*.limit.* --use-dirs'
    end
    sh 'rm -rf limits'
    sh 'mv output limits'
    sh 'touch limits/httlimit.timestamp'
  end
end

task :limit, [:chan] do |t, args|
  Rake::Task["plots/#{$jobid}/htt/#{args.chan}/limits/httlimit.timestamp"].invoke()
end

rule "httplots.timestamp" => proc{|x| x.sub('httplots', 'httlimit')} do |t|
  channel_dir = File.dirname(File.dirname(File.expand_path(t.name)))
  chan = File.basename(channel_dir)
  chdir("plots/#{$jobid}/htt/#{chan}/limits") do
    $widths.each do |width|
      sh "plotLimits.py --y-title='Coupling modifier' --x-title='M_{A} (GeV)' limits_pseudoscalar_#{width}pc.json -o limit#{width}pc --show=exp --grid --mapping='lambda x: sqrt(x)'"
    end
    sh 'touch httplots.timestamp'
  end
end

task :plotlimits, [:chan] do |t, args|
  Rake::Task["plots/#{$jobid}/htt/#{args.chan}/limits/httplots.timestamp"].invoke()
end

rule "sysbreakdown.timestamp" => proc{|x| x.sub('sysbreakdown', 'httlimit')} do |t|
  channel_dir = File.dirname(File.dirname(File.expand_path(t.name)))
  chan = File.basename(channel_dir)
  nuisances = {
    'BBB' => /_bin_[0-9]+/,
    'BTag' => /CMS_[a-z]+_b_13TeV/,
    'JESJER' => /CMS_[a-z]+_j_13TeV/,
    'QCDscale' => /QCDscale/,
    'TMass' => /TMass/,
    'Hdamp' => /Hdamp_TT/,
    'pdf' => /pdf/,
    'QCDnorm' => /QCD_[a-z]+_norm/,    
  }
  chdir("plots/#{$jobid}/htt/#{chan}/limits") do
    available=File.readlines("pseudoscalar_5pc/400/combined.txt.cmb").grep(/( lnN )|( shape )|( param )/).map {|x| x.split()[0]}
    nuisances.each do |name, group|
      to_freeze = available.select{|x| (group =~ x)}.join(',') 
      sh "combineTool.py -M Asymptotic -d */[47]*/workspace.root --there -n .limit --parallel 2 --run blind --freezeNuisances=#{to_freeze}"
      sh "combineTool.py -M CollectLimits */*/*.limit.* --use-dirs -o #{name}.json"
    end
  end  
  sh "touch #{t.name}"
end

task :sysbreakdown, [:chan] do |t, args|
  Rake::Task["plots/#{$jobid}/htt/#{args.chan}/limits/sysbreakdown.timestamp"].invoke()
  sh "python htt_scripts/sys_effects.py #{args.chan}"
end
