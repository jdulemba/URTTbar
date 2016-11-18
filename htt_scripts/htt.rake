$jobid = ENV['jobid']
$project_dir = ENV['URA_PROJECT']

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
      chdir("output") do
        sh "combineTool.py -M T2W -i pseudoscalar_#{width}pc/* -o workspace.root -P CombineHarvester.CombineTools.InterferenceModel:interferenceModel"
        sh 'combineTool.py -M Asymptotic -d */*/workspace.root --there -n .limit --parallel 5'
        sh 'combineTool.py -M CollectLimits */*/*.limit.* --use-dirs'
        sh "plotLimits.py --y-title='Coupling modifier' --x-title='M_{A} (GeV)' limits_pseudoscalar_#{width}pc.json --show=exp"
      end
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
      sh "plotLimits.py --y-title='Coupling modifier' --x-title='M_{A} (GeV)' limits_pseudoscalar_#{width}pc.json -o limit#{width}pc --show=exp"
    end
    sh 'touch httplots.timestamp'
  end
end

task :plotlimits, [:chan] do |t, args|
  Rake::Task["plots/#{$jobid}/htt/#{args.chan}/limits/httplots.timestamp"].invoke()
end
