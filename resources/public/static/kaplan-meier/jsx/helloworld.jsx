var identity = function(x) {return x;};

var get = function(key) {
  return function(o) {
    return o[key];
  }
}

var colors = [
  "rgb(200, 100, 0)",
  "rgb(100, 0, 200)",
  "rgb(0, 200, 100)"
];

var findMaxDays = function(cohort) {
  var days = cohort.individuals.map(function(individual) {return individual.days || 0})
  return Math.max.apply(null, days)
}

var findGlobalMaxDays = function(cohorts) {
  return Math.max.apply(null, cohorts.map(findMaxDays))
}

var KaplanMeierLine = React.createClass({
  render: function() {
    var line = this;

    var createFlatSegment = function(interval) {
      var xPrevious = interval.previous * line.props.dayUnits;
      var yPrevious = line.props.height - interval.remaining * line.portionUnits;
      var xDays = interval.days * line.props.dayUnits;

      return (
          <line x1={xPrevious} y1={yPrevious} x2={xDays} y2={yPrevious} stroke={interval.color} strokeWidth="3" />
      );
    }

    var createDropSegment = function(interval) {
      var xDays = interval.days * line.props.dayUnits;
      var yPrevious = line.props.height - interval.remaining * line.portionUnits;
      var yDrop = line.props.height - (interval.remaining - 1) * line.portionUnits;

      return (
          <line x1={xDays} y1={yPrevious} x2={xDays} y2={yDrop} stroke={interval.color} strokeWidth="3" />
      );
    }

    var createSegment = function(interval) {
      return (
          <g key={interval.key}>
            {createFlatSegment(interval)}
            {createDropSegment(interval)}
          </g>
      );
    }

    var expressLosses = function(cohort) {
      var losses = cohort.individuals.filter(get("days"));
      var sorted = losses.sort(function(a, b) {return a.days - b.days});
      var intervals = sorted.reduce(function(intervals, individual) {
        var lastInterval = intervals[intervals.length - 1];
        var previous = lastInterval ? lastInterval.days : 0;
        var remaining = lastInterval? lastInterval.remaining - 1 : cohort.individuals.length;
        var interval = {
          key: individual.name,
          days: individual.days,
          previous: previous,
          remaining: remaining,
          color: cohort.color
        };

        intervals.push(interval);
        return intervals;
      }, []);

      return intervals;
    }

    var cohortSize = this.props.cohort.individuals.length;
    this.portionUnits = this.props.height / cohortSize;
    var lossEvents = expressLosses(this.props.cohort);
    var lastEvent = lossEvents[lossEvents.length -1];
    var finalInterval = {
      key: "_final_",
      days: this.props.lastDay,
      previous: lastEvent.days,
      remaining: lastEvent.remaining - 1,
      color: this.props.cohort.color
    };
    var finalSegment = createFlatSegment(finalInterval);

    return (
        <g>
          {lossEvents.map(createSegment)}
          {finalSegment}
        </g>
    )
  }
});

var KaplanMeierPlot = React.createClass({
  getInitialState: function() {
    var width = 1000;
    var height = 1000;
    var tumors = ["Pancreas", "Breast", "Prostate"];
    var cohorts = {};
    var survival = [];
    
    return {
      width: width,
      height: height,
      tumors: tumors,
      cohorts: cohorts,
      survival: survival,
      dayUnits: 1,
      lastDay: 1
    };
  },

  loadSurvival: function(tumor, cohort) {
    $.ajax({
      url: this.props.survivalURL,
      dataType: 'json',
      type: 'POST',
      data: JSON.stringify(cohort),
      success: function(survival) {
        this.state.survival.push({name: tumor, individuals: survival, color: colors[this.state.survival.length]})
        var maxDays = findGlobalMaxDays(this.state.survival);
        var lastDay = Math.floor(maxDays * 1.2);
        var dayUnits = this.state.width / lastDay;

        this.setState({
          survival: this.state.survival,
          lastDay: lastDay,
          dayUnits: dayUnits
        });
      }.bind(this),
      error: function(xhr, status, err) {
        console.error(this.props.survivalURL, status, err.toString());
      }.bind(this)
    });
  },

  loadCohort: function(tumor) {
    var url = this.props.cohortURL + tumor;
    $.ajax({
      url: url,
      dataType: 'json',
      type: 'GET',
      success: function(cohort) {
        this.state.cohorts[tumor] = cohort
        this.setState({cohorts: this.state.cohorts});
        this.loadSurvival(tumor, cohort);
      }.bind(this),
      error: function(xhr, status, err) {
        console.error(url, status, err.toString());
      }.bind(this)
    });
  },

  componentDidMount: function() {
    for (var t = 0; t < this.state.tumors.length; t++) {
      this.loadCohort(this.state.tumors[t]);
    }
  },

  render: function() {
    var plot = this;

    var createLine = function(cohort) {
      return <KaplanMeierLine key={cohort.name} cohort={cohort} width={plot.state.width} height={plot.state.height} dayUnits={plot.state.dayUnits} lastDay={plot.state.lastDay} />
    }

    return (
        <div>
          <svg width={this.state.width} height={this.state.height}>
            {this.state.survival.map(createLine)}
          </svg>
        </div>
    );
  }
});

ReactDOM.render(<KaplanMeierPlot cohortURL="/gaia/individual/tumor/" survivalURL="/gaia/individual/survival" />, document.getElementById('kaplan-meier'));
