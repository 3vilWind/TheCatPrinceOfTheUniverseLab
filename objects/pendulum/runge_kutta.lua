-- RungeKutta.lua
-- ----------------------------------------------------------------- --
--      This Lua5 module is Copyright (c) 2010, Peter J Billam       --
--                        www.pjb.com.au                             --
--                                                                   --
--   This module is free software; you can redistribute it and/or    --
--          modify it under the same terms as Lua5 itself.           --
-- ----------------------------------------------------------------- --
local M = {} -- public interface
M.Version = '1.08'  -- function deepcopy is now local
M.VersionDate = '20150425'
-- 20150425 1.08 function deepcopy is now local; works with lua5.3

-- Example usage:
-- local RK = require 'RungeKutta'
-- RK.rk4()

--------------------- infrastructure ----------------------
local function arr2txt(...) -- neat printing of arrays for debug use
	local txt = {}
	for e in ... do txt[#txt+1] = string.format('%g',e) end
	return table.concat(txt,' ') .. "\n"
end
local function warn(str)
	io.stderr:write(str,'\n')
end
local function die(str)
	io.stderr:write(str,'\n')
	os.exit(1)
end
local flag = false
local a
local b
local function gaussn(standdev)
	-- returns normal distribution around 0.0 by the Box-Muller rules
	if not flag then
		a = math.sqrt(-2.0 * math.log(math.random()))
		b = 6.28318531 * math.random()
		flag = true
		return standdev * a * math.sin(b)
	else
		flag = false
		return standdev * a * math.cos(b)
	end
end
-------------------------------------------------------

function M.rk2(yn, dydt, t, dt)
	if type(yn) ~= 'table' then
		warn("RungeKutta.rk2: 1st arg must be an table\n")
		return false
	end
	if type(dydt) ~= 'function' then
		warn("RungeKutta.rk2: 2nd arg must be a function\n")
		return false
	end

	local gamma = .75;  -- Ralston's minimisation of error bounds
	local alpha = 0.5/gamma; local beta = 1.0-gamma;
	local alphadt=alpha*dt; local betadt=beta*dt; local gammadt=gamma*dt;
	local ny = #yn;
	local ynp1 = {}
	local dydtn = {}
	local ynpalpha = {}  -- Gear calls this q
	local dydtnpalpha = {}
	dydtn = dydt(t, yn);
	-- for i=1, ny do
	for i in pairs(yn) do
		ynpalpha[i] = yn[i] + alphadt*dydtn[i];
	end
	dydtnpalpha = dydt(t+alphadt, ynpalpha);
	for i in pairs(yn) do
		ynp1[i] = yn[i]+betadt*dydtn[i]+gammadt*dydtnpalpha[i];
	end
	return t+dt, ynp1
end
local function deepcopy(object)  -- http://lua-users.org/wiki/CopyTable
    local lookup_table = {}
    local function _copy(object)
        if type(object) ~= "table" then
            return object
        elseif lookup_table[object] then
            return lookup_table[object]
        end
        local new_table = {}
        lookup_table[object] = new_table
        for index, value in pairs(object) do
            new_table[_copy(index)] = _copy(value)
        end
        return setmetatable(new_table, getmetatable(object))
    end
    return _copy(object)
end


local saved_k0; local use_saved_k0 = false
function M.rk4(yn, dydt, t, dt)
	-- The Runge-Kutta-Merson 5-function-evaluation 4th-order method
	-- in the sine-cosine example, this seems to work as a 7th-order method !
	if (type(yn) ~= 'table') then
		warn("RungeKutta.rk4: 1st arg must be a table\n")
		return false
	end
	if (type(dydt) ~= 'function') then
		warn("RungeKutta.rk4: 2nd arg must be a function\n")
		return false
	end
	local ny = #yn; local i;

	local k0
	if use_saved_k0 then
		k0 = deepcopy(saved_k0)  -- a simpler single-level copy  would do...
-- without the copy() it gets trashed on the 2nd call to this function :-(
	else  k0 = dydt(t, yn)
	end
	for i in pairs(yn) do k0[i] = k0[i] * dt end

	local eta1 = {}
	for i in pairs(yn) do eta1[i] = yn[i] + k0[i]/3.0 end
	local k1 = dydt(t + dt/3.0, eta1)
	for i in pairs(yn) do k1[i] = k1[i] * dt end

	local eta2 = {}
	local k2 = {}
	for i in pairs(yn) do
		eta2[i] = yn[i] + (k0[i]+k1[i])/6.0
	end
	k2 = dydt(t + dt/3.0, eta2)
	for i in pairs(yn) do k2[i] = k2[i] * dt end

	local eta3 = {}
	for i in pairs(yn) do
		eta3[i] = yn[i] + (k0[i]+3.0*k2[i])*0.125
	end
	local k3 = dydt(t+0.5*dt, eta3)
	for i in pairs(yn) do k3[i] = k3[i] * dt end

	local eta4 = {}
	for i in pairs(yn) do
		eta4[i] = yn[i] + (k0[i]-3.0*k2[i]+4.0*k3[i])*0.5
	end
	local k4 = dydt(t+dt, eta4)
	for i in pairs(yn) do k4[i] = k4[i] * dt end

	local ynp1 = {}
	for i in pairs(yn) do
		ynp1[i] = yn[i] + (k0[i]+4.0*k3[i]+k4[i])/6.0;
	end

	-- Merson's method for error estimation, see Gear p85, only works
	-- if F is linear, ie F = Ay + bt, so that eta4 has no 4th-order
	-- errors.  So in general step-doubling is the only way to do it.
	-- Estimate error terms ...
	-- if ($epsilon) {
	-- 	my $errmax = 0; my $diff;
	-- 	for ($i=$[; $i<=$ny; $i++) {
	-- 		$diff = 0.2 * abs ($ynp1[$i] - $eta4[$i]);
	-- 		if ($errmax < $diff) { $errmax = $diff; }
	-- 	}
	-- 	-- print "errmax = $errmax\n"; -- not much related to the actual error
	-- }

	return t+dt, ynp1
end

local t = 0; local halfdt; local y2 = {}
function M.rk4_auto(yn, dydt, t, dt, arg4)
	if (type(yn) ~= 'table') then
		warn("RungeKutta.rk4_auto: 1st arg must be a table\n")
		return false
	end
	if (type(dydt) ~= 'function') then
		warn("RungeKutta.rk4_auto: 2nd arg must be a function\n")
		return false
	end
	if dt == 0 then dt = 0.1 end
	local errors; local epsilon = nil
	if (type(arg4) == 'table') then
		errors = arg4; epsilon = nil
	else
		epsilon = math.abs(arg4); errors = nil
		if epsilon == 0 then epsilon = .0000001 end
	end
	local ny = #yn; local i

	local y1 = {}
	local y3 = {}
	saved_k0 = dydt(t, yn)
	local resizings = 0;
	local highest_low_error = 0.1e-99; local highest_low_dt = 0.0;
	local lowest_high_error = 9.9e99;  local lowest_high_dt = 9.9e99;
	while true do
		halfdt = 0.5 * dt; local dummy
		use_saved_k0 = true
		dummy, y1 = M.rk4(yn, dydt, t, dt)
		dummy, y2 = M.rk4(yn, dydt, t, halfdt)
		use_saved_k0 = false
		dummy, y3 = M.rk4(y2, dydt, t+halfdt, halfdt)

		local relative_error
		if epsilon then
	 		local errmax = 0; local diff; local ymax = 0
	 		for i in pairs(yn) do
	 			diff = math.abs(y1[i] - y3[i])
	 			if errmax < diff then errmax = diff end
	 			if ymax < math.abs(yn[i]) then ymax = math.abs(yn[i]) end
	 		end
			relative_error = errmax / (epsilon*ymax)
		elseif errors then
			relative_error = 0.0; local diff;
	 		for i in pairs(yn) do
	 			diff = math.abs(y1[i] - y3[i]) / math.abs(errors[i])
	 			if relative_error < diff then relative_error = diff end
	 		end
		else
			die "RungeKutta.rk4_auto: epsilon & errors both undefined\n";
		end
		-- Gear's "correction" assumes error is always in 5th-order terms :-(
		-- $y1[$i] = (16.0*$y3{$i] - $y1[$i]) / 15.0;
		if relative_error < 0.60 then
			if dt > highest_low_dt then
				highest_low_error = relative_error; highest_low_dt = dt
			end
		elseif relative_error > 1.67 then
			if dt < lowest_high_dt then
				lowest_high_error = relative_error; lowest_high_dt = dt
			end
		else
			break
		end
		if lowest_high_dt<9.8e99 and highest_low_dt>1.0e-99 then -- interpolate
			local denom = math.log(lowest_high_error/highest_low_error)
			if highest_low_dt==0.0 or highest_low_error==0.0 or denom == 0.0 then
				dt = 0.5 * (highest_low_dt+lowest_high_dt)
			else
				dt = highest_low_dt * ( (lowest_high_dt/highest_low_dt)
				 ^ ((math.log(1.0/highest_low_error)) / denom) )
			end
		else
			local adjust = relative_error^(-0.2) -- hope error is 5th-order ...
			if math.abs(adjust) > 2.0 then
				dt = dt * 2.0  -- prevent infinity if 4th-order is exact ...
			else
				dt = dt * adjust
			end
		end
		resizings = resizings + 1
		if resizings>4 and highest_low_dt>1.0e-99 then
			-- hope a small step forward gets us out of this mess ...
			dt = highest_low_dt;  halfdt = 0.5 * dt;
			use_saved_k0 = true
			dummy, y2 = M.rk4(yn, dydt, t, halfdt)
			use_saved_k0 = false
			dummy, y3 = M.rk4(y2, dydt, t+halfdt, halfdt)
			break
		end
	end
	return t+dt, dt, y3
end

function M.rk4_auto_midpoint()
	return t+halfdt, y2
end

------------------------ EXPORT_OK routines ----------------------

function M.rk4_ralston (yn, dydt, t, dt)
	if (type(yn) ~= 'table') then
		warn("RungeKutta.rk4_ralston: 1st arg must be arrayref\n")
		return false
	end
	if (type(dydt) ~= 'function') then
		warn("RungeKutta.rk4_ralston: 2nd arg must be a subroutine ref\n")
		return false
	end
	local ny = #yn; local i;

	-- Ralston's minimisation of error bounds, see Gear p36
	local alpha1=0.4; local alpha2 = 0.4557372542  -- = .875 - .1875*(sqrt 5);

	local k0 = dydt(t, yn)
	for i in pairs(yn) do k0[i] = dt * k0[i] end

	local k1 = {}
	for i in pairs(yn) do k1[i] = yn[i] + 0.4*k0[i] end
	k1 = dydt(t + alpha1*dt, k1)
	for i in pairs(yn) do k1[i] = dt * k1[i] end

	local k2 = {}
	for i in pairs(yn) do
		k2[i] = yn[i] + 0.2969776*k0[i] + 0.15875966*k1[i]
	end
	k2 = dydt(t + alpha2*dt, k2)
	for i in pairs(yn) do k2[i] = dt * k2[i] end

	local k3 = {}
	for i in pairs(yn) do
		k3[i] = yn[i] + 0.21810038*k0[i] - 3.0509647*k1[i] + 3.83286432*k2[i]
	end
	k3 = dydt(t+dt, k3)
	for i in pairs(yn) do k3[i] = dt * k3[i] end

	local ynp1 = {}
	for i in pairs(yn) do
		ynp1[i] = yn[i] + 0.17476028*k0[i]
		 - 0.55148053*k1[i] + 1.20553547*k2[i] + 0.17118478*k3[i]
	end
	return t+dt, ynp1
end
function M.rk4_classical(yn, dydt, t, dt)
	if (type(yn) ~= 'table') then
		warn("RungeKutta.rk4_classical: 1st arg must be arrayref\n")
		return false
	end
	if (type(dydt) ~= 'function') then
		warn("RungeKutta.rk4_classical: 2nd arg must be subroutine ref\n")
		return false
	end
	local ny = #yn; local i;

	-- The Classical 4th-order Runge-Kutta Method, see Gear p35

	local k0 = dydt(t, yn)
	for i in pairs(yn) do k0[i] = dt * k0[i] end

	local eta1 = {}
	for i in pairs(yn) do eta1[i] = yn[i] + 0.5*k0[i] end
	local k1 = dydt(t+0.5*dt, eta1)
	for i in pairs(yn) do k1[i] = dt * k1[i] end

	local eta2 = {}
	for i in pairs(yn) do eta2[i] = yn[i] + 0.5*k1[i] end
	local k2 = dydt(t+0.5*dt, eta2)
	for i in pairs(yn) do k2[i] = dt * k2[i] end

	local eta3 = {}
	for i in pairs(yn) do eta3[i] = yn[i] + k2[i] end
	local k3 = dydt(t+dt, eta3)
	for i in pairs(yn) do k3[i] = dt * k3[i] end

	local ynp1 = {}
	for i in pairs(yn) do
		ynp1[i] = yn[i] + (k0[i] + 2.0*k1[i] + 2.0*k2[i] + k3[i]) / 6.0;
	end
	return t+dt, ynp1
end


return M
