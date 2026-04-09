import numpy as np

file_path = "msd_results.txt"  

def read_msd_file(file_path):
    steps = []
    msd_values = []
    
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("Step:"):
                parts = line.split()
                steps.append(int(parts[1]))  
                msd_values.append(float(parts[9]))  
    return steps, msd_values


def calculate_diffusion_direct(steps, msd_values, timestep):
    
    time_intervals = [(steps[i+1] - steps[i]) * timestep for i in range(len(steps) - 1)]
    msd_increments = [msd_values[i+1] - msd_values[i] for i in range(len(msd_values) - 1)]
    
    avg_msd_increment = sum(msd_increments) / len(msd_increments)
    avg_time_interval = sum(time_intervals) / len(time_intervals)
    
    diffusion_coefficient = avg_msd_increment / (4 * avg_time_interval)
    return diffusion_coefficient, avg_msd_increment, avg_time_interval


def main():
    timestep = 0.01
    steps, msd_values = read_msd_file(file_path)

    diffusion_coefficient, avg_msd_increment, avg_time_interval = calculate_diffusion_direct(steps, msd_values, timestep)

    print(f"diffusion_coefficient: {diffusion_coefficient:.5f}")
    print(f"avg_msd: {avg_msd_increment:.5f}")

if __name__ == "__main__":
    main()
    
