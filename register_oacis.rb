repo_dir = File.expand_path(File.dirname(__FILE__))

localhost = Host.find_by_name("localhost")

sim_params = {
  name: "NetworkSamplingFixedF0",
  command: "ruby #{repo_dir}/runner/run_power_mean_sampling.rb _input.json",
  support_input_json: true,
  print_version_command: "cd #{repo_dir} && git describe --always",
  parameter_definitions: [
    {key: "N", type: "Integer", default: 1000, description: "network size"},
    {key: "k", type: "Float", default: 20.0, description: "original degree"},
    {key: "f0", type: "Float", default: 0.3, description: "typical scale of f"},
    {key: "alpha", type: "Float", default: 1.0, description: "exponent of the Weibull distribution"},
    {key: "beta", type: "Float", default: 0.0, description: "exponent of the generalized mean"}
  ],
  description: "Network sampling model with fixed f0",
  executable_on: [ localhost ]
}

sim_params2 = {
  name: "NetworkSamplingTunedF0",
  command: "ruby #{repo_dir}/runner/run_power_mean_sampling_tuned_f0.rb _input.json",
  support_input_json: true,
  print_version_command: "cd #{repo_dir} && git describe --always",
  parameter_definitions: [
    {key: "N", type: "Integer", default: 1000, description: "network size"},
    {key: "k0", type: "Float", default: 150.0, description: "original degree"},
    {key: "expected_k", type: "Float", default: 15.0, description: "sampled degree"},
    {key: "dk", type: "Float", default: 2.0, description: "allowed difference of the sampled degree"},
    {key: "alpha", type: "Float", default: 1.0, description: "exponent of the Weibull distribution"},
    {key: "beta", type: "Float", default: 0.0, description: "exponent of the generalized mean"}
  ],
  description: "Network sampling model with tuned f0",
  executable_on: [ localhost ]
}

analyzer_params = {
  name: "make_plot",
  command: "#{repo_dir}/network_analysis/plot/plot_all.sh _input .",
  support_input_json: true,
  type: "on_run",
  auto_run: "first_run_only",
  print_version_command: "cd #{repo_dir} && git describe --always",
  description: "make plot of network properties.",
  executable_on: [ localhost ],
  auto_run_submitted_to: localhost
}

if Simulator.where(name: sim_params[:name]).exists?
	puts "simulator #{sim_params[:name]} already exists" 
else
	sim = Simulator.create!(sim_params)
  sim.analyzers.create!(analyzer_params)
end

if Simulator.where(name: sim_params2[:name]).exists?
	puts "simulator #{sim_params2[:name]} already exists" 
else
	sim2 = Simulator.create!(sim_params2)
  sim2.analyzers.create!(analyzer_params)
end

