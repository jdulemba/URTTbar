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
