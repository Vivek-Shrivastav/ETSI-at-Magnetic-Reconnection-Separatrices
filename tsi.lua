-- Gkyl ------------------------------------------------------------------------
-- 1x1v simulation of two stream instability
--------------------------------------------------------------------------------
-- App dependencies
--------------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()
--------------------------------------------------------------------------------
-- Preamble
--------------------------------------------------------------------------------
local vth_e = 0.3             -- electron thermal velocity
local drift = 1.0             -- drift velocity
local perturbK = 0.5          -- perturbation wave-number
local perturbMag = 1.0e-3     -- distribution function perturbation
local function maxwellian(v, n, u, vth)
   return n/(math.sqrt(2*math.pi)*vth)
      *math.exp(-0.5*((v-u)/vth)^2)
end
plasmaApp = Plasma.App {
   -----------------------------------------------------------------------------
   -- Common
   -----------------------------------------------------------------------------
   logToFile = true,
   tEnd = 1000.0,                  -- end time
   nFrame = 1000,                 -- number of output frames
   lower = {-2*math.pi/perturbK},  -- configuration space lower left
   upper = {2*math.pi/perturbK},   -- configuration space upper right
   cells = {512},                 -- number of configuration space cells
   basis = "serendipity",        -- basis function
   polyOrder = 2,                -- polynomial order
   timeStepper = "rk3",          -- 3rd order Runge-Kutta
   -- MPI decomposition for configuration space
   decompCuts = {64},   -- number of cuts in each configuration direction
   useShared = false,  -- flag for MPI shared memory
   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions
   -----------------------------------------------------------------------------
   -- Electrons
   -----------------------------------------------------------------------------
   elc = Plasma.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-10.0},
      upper = {10.0},
      cells = {1024},
      decompCuts = {1},
      -- initial conditions
      init = function (t, z)
	 local x, vx = z[1], z[2]
	 local f = maxwellian(vx, 0.5, drift, vth_e) +
            maxwellian(vx, 0.5, -drift, vth_e)
	 return (1+perturbMag*math.cos(perturbK*x))*f
      end,
      evolve = true,
      diagnostics = {"M0","M1i","M2ij","M3i","intM0","intM1i","intM2Flow","intM2Thermal"},
   },
   
   -----------------------------------------------------------------------------
   -- Field Solver
   -----------------------------------------------------------------------------
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, z)
         local x = z[1]
	 local Ex = -perturbMag*math.sin(perturbK*x)/perturbK
	 return Ex, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true,
   },
}
--------------------------------------------------------------------------------
-- Run the application
--------------------------------------------------------------------------------
plasmaApp:run()


