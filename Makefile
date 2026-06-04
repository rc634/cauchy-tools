CXX      = g++
CXXFLAGS = -std=c++17 -O2 -MMD -MP
TARGET   = ahfinder
BUILDDIR = build

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:%.cpp=$(BUILDDIR)/%.o)
DEPS := $(SRCS:%.cpp=$(BUILDDIR)/%.d)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILDDIR)/%.o: %.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

-include $(DEPS)

.PHONY: clean
clean:
	rm -rf $(BUILDDIR) $(TARGET)
