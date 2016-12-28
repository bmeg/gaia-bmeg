function Pie() {
  var slice = d3.layout.pie().value(function(d) {return d.population});
  var arc = d3.svg.arc().innerRadius(30).outerRadius(200);
  var color = d3.scale.category10();

  var svg = d3.select('svg.cohort');
  var p = svg.append('g').attr('transform', 'translate(400, 250)');
  var pie = p.selectAll('path.slice');
  var l = svg.append('g').attr('class', 'legend')
  var legend = l.selectAll('text.line');

  function pieLayout(cohort, value) {
    var total = cohort.reduce(function(t, d) {return t + value(d)}, 0);
    return cohort.reduce(function(l, d) {
      var v = value(d);
      var start = l.length > 0 ? l[l.length - 1].endAngle : 0;
      var slice = {
        startAngle: start,
        endAngle: start + v * 2 * Math.PI / total,
        padAngle: 0,
        value: v,
        data: d
      };

      l.push(slice);
      return l;
    }, []);
  }

  function update(cohort) {
    cohort = cohort.sort(function(a, b) {return (a.tumor < b.tumor) ? -1 : (a.tumor > b.tumor) ? 1 : 0});
    var slices = pieLayout(cohort, function(d) {return d.population});

    pie = pie.data(slices, function(d, i) {return d.data.tumor});
    pie.exit().remove();
    pie.enter()
      .append('path')
      .attr('class', 'slice')
      .attr('fill', function(d) {return color(d.data.tumor)});

    pie.transition()
      .duration(1000)
      .attr('d', arc)

    legend = legend.data(slices, function(d, i) {return d.data.tumor});
    legend.exit().remove();
    legend.enter()
      .append('text')
      .attr('class', 'line')
      .attr('fill', function(d) { return color(d.data.tumor); })
      .text(function(d) { return '• ' + d.data.tumor; })

    legend.transition()
      .duration(1000)
      .attr('y', function(d, i) { return 20 * (i + 1); });
  }

  update([]);

  return {
    update: update
  }
}

var GeneInput = React.createClass({
  getInitialState: function() {
    return {
      input: "",
      lastMatch: "none",
      tumorCounts: {},
      pie: Pie()
    };
  },

  changeInput: function(event) {
    var gene = event.target.value
    this.setState({input: gene})

    var url = "/gaia/gene/" + gene + "/tumor/counts";

    $.ajax({
      url: url,
      dataType: 'json',
      type: 'GET',
      success: function(result) {
        if (Object.keys(result).length > 0) {
          this.setState({tumorCounts: result, lastMatch: this.state.input})
          var slices = Object.keys(result).map(function(key) {
            return {'tumor': key, 'population': result[key]};
          });

          this.state.pie.update(slices)
        }
      }.bind(this),

      error: function(xhr, status, err) {
        console.error(url, status, err.toString());
      }.bind(this)
    });
  },

  render: function() {
    return (
      <div>
        <form onChange={(e) => this.changeInput(e)}>
          <div className="mdl-textfield mdl-js-textfield mdl-textfield--floating-label">
            <input id="gene-text-input" type="text" name="gene" className="mdl-textfield__input" />
            <label className="mdl-textfield__label" htmlFor="gene-text-input">Gene...</label>
          </div>
        </form>
        <p>Showing tumor sites for variants in: {this.state.lastMatch}</p>
      </div>
    );
  }
});

ReactDOM.render(<GeneInput />, document.getElementById('gene-input'));
