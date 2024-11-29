#include <future>
#include <iomanip>
#include <iostream>
#include <random>
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include <__random/random_device.h> /* Hardware random number generator */
#include <semaphore>

// Forward declaration to allow the RunSimulation function to call this function
std::pair<int, int> RunSubsetSimulation(int numSimulations, int roomsPerLevelLayout[], int numLayouts, int coinsPerRoom[], int numRooms);

// Forward declaration to allow the GenerateLevelAndGetCoins function to call this function
int GenerateLevelAndGetCoins(int levelLayoutType, int roomsPerLevelLayout[], int coinsPerRoom[], int numRooms);

// Factorial: N! = N * (N-1) * (N-2) * ... * 1
// Method: Multiply all the numbers from N down to 1
// 5! = 5 * 4 * 3 * 2 * 1 = 120
double Factorial(double n) {
    // Convert n to an integer for the loop
    const auto number = static_cast<int>(n);

    if (number <= 1)
        return 1.0;

    // Start with 1
    double result = 1.0;

    // Multiply by each number from n down to 1
    for (int i = number; i > 1; i--)
    {
        // Multiply the result by the current number
        result *= i;
    }

    //std::cout << "Factorial of " << n << " is " << std::fixed << std::setprecision(0) << result << "\n";

    return result;
}

// CalculateNumPosibleLayouts: Calculate the number of possible layouts
// Method: Calculate the number of possible layouts using permutations or combinations
// Permutations: The order matters, so the number of possible layouts is P(objects, sample) = objects! / (objects - sample)!
// Combinations: The order does not matter, so the number of possible layouts is C(objects, sample) = objects! / (sample! * (objects - sample)!)
// 5 rooms, 3 selected: C(5, 3) = 5! / (3! * (5 - 3)!) = 10
// 5 rooms, 3 selected: P(5, 3) = 5! / (5 - 3)! = 60
// If the order does not matter, use combinations, else use permutations
double CalculateNumPosibleLayouts(double sample, double objects, const bool useCombinations = false)
{
    // Check for invalid cases
    if (sample > objects || sample < 0 || objects < 0)
        return 0.0; // Invalid case

    // Check if we are using permutations
    const bool usePermutations = !useCombinations;

    if (useCombinations)
    {
        return Factorial(objects) / (Factorial(sample) * Factorial(objects - sample));
    }

    if (usePermutations)
    {
        return Factorial(objects) / Factorial(objects - sample);
    }
}

// RunSimulation: Run the Monte Carlo simulation
// Method: Run a subset of the total number of simulations in parallel using multiple threads
// Calculate the total coins and the highest coins from the subset simulations
// Accumulate the results from the subset simulations
// Calculate the average, ratio, and duration of the simulation
void RunSimulation(int numSimulations, int roomsPerLevelLayout[], int numLayouts, int coinsPerRoom[], int numRooms)
{
    // Get the number of concurrent threads supported by the hardware
    unsigned int n = std::thread::hardware_concurrency();
    std::cout << "\n" << n << " Concurrent threads are supported.\n";

    // Record the start time
    auto startTime = std::chrono::high_resolution_clock::now();

    // Limit the number of simulations per thread to avoid running out of memory
    // Constexpr is used to define a constant expression at compile time
    constexpr int maxSimulationsPerThread = 100000;

    // Calculate the number of threads needed (Ceiling division)
    int threads = (numSimulations + maxSimulationsPerThread - 1) / maxSimulationsPerThread;

    // Calculate the number of simulations per thread
    int simulationsPerThread = numSimulations / threads;

    // Calculate the remainder of simulations
    int remainder = numSimulations % threads;

    // Start a new thread and store the future results in the vector
    std::vector<std::future<std::pair<int, int>>> futures;

    for (int i = 0; i < threads; ++i)
    {
        // Calculate the number of simulations for this thread
        int currentSimulations = simulationsPerThread + (i < remainder ? 1 : 0);

        futures.push_back(
                        std::async(
                                    std::launch::async,         // Run the function in a separate thread
                                    RunSubsetSimulation,
                                    currentSimulations,
                                    roomsPerLevelLayout,
                                    numLayouts,
                                    coinsPerRoom,
                                    numRooms));
    }

    // Accumulate the results from the futures
    long totalCoins = 0;    // Use long to avoid overflow while testing 100,000,000 simulations
    int highestCoins = 0;

    // Wait for all threads to finish and get the results
    for (auto& future : futures)
    {
        // Get the result from the future
        // Using auto here because we can create and assign a pair of integers in one line
        auto [subsetTotalCoins, subsetMaxCoins] = future.get();

        // Add the subset results to the total
        totalCoins += subsetTotalCoins;

        // Check if the subset max coins is higher than the current highest coins
        // If it is, update the highest coins
        if (subsetMaxCoins > highestCoins)
            highestCoins = subsetMaxCoins;
    }

    // Calculate the average
    double averageCoins = static_cast<double>(totalCoins) / numSimulations;

    // Calculate the ratio
    double ratio = averageCoins / highestCoins;

    // Record the end time
    auto endTime = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    std::chrono::duration<double> elapsedSeconds = endTime - startTime;

    std::cout << "\nMonte Carlo Simulation Results:\n";
    std::cout << "Simulation duration: " << elapsedSeconds.count() << " seconds\n";
    std::cout << "Total number of simulations: " << numSimulations << "\n";
    std::cout << "Total number of coins: " << totalCoins << "\n";
    std::cout << "Highest number of coins: " << highestCoins << "\n";
    std::cout << "Average number of coins: " << averageCoins << "\n";
    std::cout << "Ratio (decimal): " << ratio << "\n";
    std::cout << "Ratio (percentage): " << std::fixed << std::setprecision(1) << ratio * 100 << "%\n";
    std::cout << "Ratio (fraction, tenths): " << static_cast<int>(std::round(ratio * 10)) << "/10\n";}

// RunSubsetSimulation: Run a subset of the total number of simulations
// Method: Randomly choose a layout index
// Generate the level and get the coins
// Add the coins to the total
// Check if the coins are higher than the current highest coins
// If they are, update the highest coins
std::pair<int, int> RunSubsetSimulation(int numSimulations, int roomsPerLevelLayout[], int numLayouts, int coinsPerRoom[], int numRooms)
{
    // Random number generation setup
    // std::random_device is a hardware random number generator
    std::random_device randDevice;

    // std::mt19937 is a Mersenne Twister pseudo-random generator
    // Seeded with the hardware random number generator
    std::mt19937 random(randDevice());

    // std::uniform_int_distribution is a random number distribution that produces integers with equal probability
    std::uniform_int_distribution<int> layoutDist(0, numLayouts - 1);

    // Accumulate the results from the simulations
    int totalCoins = 0;
    int highestCoins = 0;

    for (int i = 0; i < numSimulations; ++i)
    {
        // Randomly choose a layout index
        const int layoutIndex = layoutDist(random);

        // Generate the level and get the coins
        const int coins = GenerateLevelAndGetCoins(layoutIndex, roomsPerLevelLayout, coinsPerRoom, numRooms);

        // Add the coins to the total
        totalCoins += coins;

        // Check if the coins are higher than the current highest coins
        // If they are, update the highest coins
        if (coins > highestCoins)
            highestCoins = coins;
    }

    // Return the int pair of total coins and highest coins
    return {totalCoins, highestCoins};
}

int GenerateLevelAndGetCoins(int levelLayoutType, int roomsPerLevelLayout[], int coinsPerRoom[], int numRooms)
{
    // Track wether a room was used or not, as we can only use each room once
    bool* roomsUsed = new bool[numRooms];
    for (int i = 0; i < numRooms; ++i) roomsUsed[i] = false; // Set all rooms to NOT used by default

    // Track coins total in all rooms, start at 0
    int coinsTotal = 0;

    // For each room in the chosen layout...
    for (int i = 0; i < roomsPerLevelLayout[levelLayoutType]; ++i)
    {
        // Start with 0 room selected (will be overridden)
        int roomSelected = 0;
        do
        {
            // Randomly choose a room
            roomSelected = rand() % numRooms;
        } while (roomsUsed[roomSelected] == true); // Keep choosing a room until you get one that hasn't been used

        // Mark the room as used
        roomsUsed[roomSelected] = true;

        // Add the coins for the randomly selected room to the total coins
        coinsTotal += coinsPerRoom[roomSelected];
    }

    // Return the total coins for this room layout
    return coinsTotal;
}

void CalculateNumLayoutsForAllTypes(int roomsPerLevelLayout[], int numLayouts, int numRooms, bool verbose = false)
{
    std::cout << "Number of layouts possible, by layout type:\n";
    std::cout << "Showing results in " << (verbose ? "fixed" : "scientific") <<" notation.\n";

    // Loop through all layout types
    for (int i = 0; i < numLayouts; ++i) {
        // Add line breaks every 5 layout types
        if (i % 5 == 0)
        {
            std::cout << "\n |\t";
        }
        // Calculate the possible layouts for this number of rooms
        double numLayouts = CalculateNumPosibleLayouts(roomsPerLevelLayout[i], numRooms);

        // By default, output in scientific notation
        // For testing I wanted to see the full number

        if (verbose)
            std::cout << "Layout " << (i+1) << ": " << std::fixed << std::setprecision(0) << numLayouts << "\t|\t";
        else
            std::cout << "Layout " << (i+1) << ": " << numLayouts << "\t|\t"; // Scientific notation
    }
    std::cout << "\n";

}

int main()
{
    // random number generation setup
    srand(time(0));

    // Room layouts
    int roomsPerLevelLayoutType[] =
    { 10, 10, 10, 10, 11, 11, 12, 12, 12, 13,
    14, 14, 14, 14, 14, 14, 14, 15, 15, 15,
    15, 16, 16, 16, 16, 17, 17, 17, 18, 18,
    18, 18, 19, 19, 20, 20, 20, 20, 20, 21,
    22, 22, 22, 23, 23, 24, 24, 26, 27, 30
    };
    int numLayoutTypes = 50;

    // Coins
    int coinsPerRoom[] =
    {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
        7, 7, 7, 7, 8, 8, 8, 9, 9, 10
    };
    int numRooms = 120;
    int numSimulations = 10000000; // 10,000,000

    CalculateNumLayoutsForAllTypes(roomsPerLevelLayoutType, numLayoutTypes, numRooms);
    std::cout << "\nRunning simulation, " << numSimulations << " iterations, please wait...";

    RunSimulation(numSimulations, roomsPerLevelLayoutType, numLayoutTypes, coinsPerRoom, numRooms);
}
