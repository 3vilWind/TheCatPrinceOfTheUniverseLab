local utils = require "objects.pendulum.utils"

local M = {}

local objects = {}

function M.get_sample(ind, sample_num)
	local obj = objects[ind]
	return obj:get_sample(sample_num)
end

function M.get_index(obj)
	for k, v in pairs(objects) do
		if v:equals(obj) then
			return k
		end
	end
	local ind = math.random(1, 1000000)
	objects[ind] = obj
	return ind
end


return M
