local utils = require "objects.pendulum.utils"
local RK = require "objects.pendulum.runge_kutta"

local SimplePendulum = {}

function SimplePendulum:new(sample_interval, initial_angle, period, g)
	local obj = {}
	setmetatable(obj, self)
	self.__index = self
	obj.samples = {}
	obj.sample_interval = sample_interval
	obj.type = 1
	obj.initial_angle = initial_angle
	obj.period = period
	obj.g = g
	local m = utils.agm(1, math.cos(initial_angle/2))
	obj.rope_length = (period * period * m * m * g) / (4 * math.pi * math.pi)
	obj.last_sample = 0
	obj.last_t = 0
	obj.last_y = {obj.initial_angle, 0}

	function f(t, y)
		local theta = y[1]
		local x = y[2]
		return {x, -(obj.g/obj.rope_length) * math.sin(theta)}
	end
	
	obj.f = f
	return obj
end

function SimplePendulum:equals(other)
	return self.type == other.type and
	self.sample_interval == other.sample_interval and
	self.initial_angle == other.initial_angle and
	self.period == other.period and
	self.g == other.g
end

function SimplePendulum:get_sample(sample_num)
	for i = self.last_sample, sample_num do
		self:compute_next_sample()
	end

	return self.samples[sample_num]
end

function SimplePendulum:compute_next_sample()
	self.last_t, self.last_y = RK.rk4(self.last_y, self.f, self.last_t, self.sample_interval)
	self.samples[self.last_sample] = self.last_y[1]
	self.last_sample = self.last_sample + 1
end

return SimplePendulum
