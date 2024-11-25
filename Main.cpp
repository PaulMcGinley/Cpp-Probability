#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

double Factorial(double n)
{
    // INSERT CODE HERE
    // Write this helper function to assist in writing the CalculateNumPossibleLayouts function below
    return 0.0;
}

double CalculateNumPosibleLayouts(double sample, double objects)
{
    // INSERT CODE HERE
    // Determine whether the randomly generated elements will be combinations or permutations
    // Based on that information, write code to calculate the number of unique layouts possible on each level type
    // This function only needs to calculate for a single layout type, it will be called by a larger function below.
    // You should write and use the Factorial funciton above.
    return 0.0;
}


int GenerateLevelAndGetCoins(int levelLayoutType, int roomsPerLevelLayout[], int coinsPerRoom[], int numRooms);
void RunSimulation(int numSimulations, int roomsPerLevelLayout[], int numLayouts, int coinsPerRoom[], int numRooms)
{
    // INSERT CODE HERE
    // Write code to run a Monte Carlo simulation over 10,000,000 iterations to determine the highest number of 
    // coins that can reasonably be generated, as well as the average number of coins.
    // Finally, calculate the ratio between the average number of coins, and the expected maximum number of coins, 
    // and display it in three forms: a decimal ratio, a percentage to one decimal place, and as a fraction (tenths).
    // You may use the GenerateLevelAndGetCoins function above to assist in your simulation - be sure
    // to randomly choose a layout index first and pass it to the function!
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

void CalculateNumLayoutsForAllTypes(int roomsPerLevelLayout[], int numLayouts, int numRooms)
{
    std::cout << "Number of layouts possible, by layout type:";
    // Loop through all layout types
    for (int i = 0; i < numLayouts; ++i)
    {
        // Add line breaks every 5 layout types
        if (i % 5 == 0)
        {
            std::cout << "\n |\t";
        }
        // Calculate the possible layouts for this number of rooms
        double numLayouts = CalculateNumPosibleLayouts(roomsPerLevelLayout[i], numRooms);
        std::cout << "Layout " << (i+1) << ": " << numLayouts << "\t|\t";
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
    std::cout << "\n\nRunning simulation, " << numSimulations << "iterations, please wait...";
    RunSimulation(numSimulations, roomsPerLevelLayoutType, numLayoutTypes, coinsPerRoom, numRooms);
}
