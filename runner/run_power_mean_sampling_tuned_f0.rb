require 'json'
require 'date'
require 'fileutils'
require 'logger'

# Usage:
#   ruby #{__FILE__} _input.json

# execute simulator
unless ARGV[0]
  $stderr.puts "Usage: ruby #{File.basename(__FILE__)} _input.json"
  raise "invalid argument"
end

# calculate f_0 so that we have expected degree in the sampled network

params = JSON.load(File.open(ARGV[0]))

simulator = File.expand_path( File.join( File.dirname(__FILE__), "../simulator/power_mean_sampling_truncation_tuned_f0.out") )

$stderr.puts "Running simulation"
keys = %w(N k0 f0 expected_k dk alpha power _seed)
args = keys.map {|key| params[key] }
command = "#{simulator} #{args.join(' ')}"
$stderr.puts "Running simulator : #{DateTime.now}"
$stderr.puts command
system(command)
raise "Simulator failed" unless $?.to_i == 0

sleep 1

# execute analyzer
analyzer = File.expand_path( File.join( File.dirname(__FILE__), "../analyzer/analyzer.out") )
edge_file = "sampled.edg"
command = "#{analyzer} #{edge_file}"
$stderr.puts "Running analyzer : #{DateTime.now}"
$stderr.puts command
system(command)
raise "Analyzer failed" unless $?.to_i == 0

# load f0
f0 = File.read("f0.dat").to_f
json = JSON.load( File.open("_output.json") )
json["f0"] = f0
File.open("_output.json", 'w') {|io| io.puts json.to_json }


