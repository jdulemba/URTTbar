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

ttxs_rules = "#{$project_dir}/ttxs_scripts/ttxs.rake"
import ttxs_rules

ctag_rules = "#{$project_dir}/ctag_scripts/ctag.rake"
import ctag_rules

htt_rules = "#{$project_dir}/htt_scripts/htt.rake"
import htt_rules
