# Message for Khadija 💕

## Running the Twisted Bilayer Graphene Band Structure Calculator

Hey beautiful! I've converted those Fortran files to Python for your professor's research project. Here's everything you need to know to run the band structure calculations:

### 📋 **Quick Setup**

1. **Install the required packages:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Test that everything works:**
   ```bash
   python test_conversion.py
   ```
   You should see "All tests passed!" at the end.

3. **Run the main calculation:**
   ```bash
   python plot_bands.py
   ```

### 🔬 **What This Does**

This code calculates the electronic band structure for **twisted bilayer graphene** - a fascinating 2D material where two graphene sheets are rotated relative to each other. The physics becomes really interesting at "magic angles" (~1.1°) where superconductivity can emerge!

The calculation:
- Sets up a moiré superlattice for the twisted layers
- Builds the tight-binding Hamiltonian with inter-layer coupling
- Applies lattice relaxation (atomic positions adjust due to interlayer forces)
- Calculates bands along the Γ-K-M-Γ path in the Brillouin zone
- Outputs energy eigenvalues relative to the neutrality point

### ⚙️ **Key Parameters** (in `hamiltonian.py`)

- **`NTHETA = 5`** - Controls the twist angle (currently ~1.05°)
- **`NUMK = 30`** - Number of k-points for sampling (higher = more detailed but slower)
- **`NRELAX = 1`** - Enables lattice relaxation (set to 0 to disable)
- **`PRESSURE = 1.0`** - Hydrostatic pressure parameter

### 📊 **Output**

The program creates a file like `Bands-I5-numk30-relax1.dat` containing:
- Column 1: k-point coordinate along the path (0 to 1)
- Columns 2+: Band energies in eV (relative to neutrality point)

### 🎯 **For Your Professor**

This Python version should give **identical results** to the original Fortran code while being:
- ✅ Easier to modify and extend
- ✅ Cross-platform compatible
- ✅ Well-documented with type hints
- ✅ Parallelized using modern Python libraries

The conversion preserves all the sophisticated physics including:
- Complete tight-binding model with realistic parameters
- Lattice relaxation coefficients for multiple twist angles
- Proper treatment of the moiré superlattice geometry
- Efficient partial matrix diagonalization

### 🚀 **Advanced Usage**

To modify the twist angle, change `NTHETA` in `hamiltonian.py`. Note that relaxation coefficients are only available for `NTHETA = 2-35`. For new angles, you'd need to compute new relaxation parameters.

### 💡 **Troubleshooting**

- **Import errors**: Make sure all dependencies are installed with `pip install -r requirements.txt`
- **Slow performance**: The code uses all CPU cores by default. For large systems, consider reducing `NUMK` for testing
- **Memory issues**: Reduce `NTHETA` or `NUMK` if running on limited memory systems

Hope this helps with your research! The physics here is absolutely beautiful - twisted bilayer graphene is one of the most exciting discoveries in condensed matter physics in recent years. 

Let me know if you or your professor have any questions about the conversion or need help extending the code!

xoxo ❤️ 