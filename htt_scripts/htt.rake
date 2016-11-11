$jobid = ENV['jobid']

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

task :limit,[:chan] do |t, args|
  chdir("plots/#{$jobid}/htt/#{args.chan}") do
    sh "setup_andrey.py #{args.chan}.root"
    sh 'rm -rf limits'
    sh 'mv output limits'
    chdir("limits") do
      sh 'combineTool.py -M T2W -i pseudoscalar/* -o workspace.root -P CombineHarvester.CombineTools.InterferenceModel:interferenceModel'
      sh 'combineTool.py -M Asymptotic -d */*/workspace.root --there -n .limit --parallel 4'
      sh 'combineTool.py -M CollectLimits */*/*.limit.* --use-dirs -o limits.json'
      sh 'plotLimits.py --y-title="Coupling modifier" --x-title="M_{A} (GeV)" limits_pseudoscalar.json --show=exp'
    end
  end
end
