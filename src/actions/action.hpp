#ifndef __ACTION_HPP__
#define __ACTION_HPP__

class Action{
	public:
		virtual void operator()();
		virtual void configure();
		virtual ~Action();
};

#endif /* __ACTION_HPP__ */
