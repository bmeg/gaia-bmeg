var snipPrefix = function(s) {
  return s.substring(s.indexOf(':') + 1);
}

var exploreVertex = function(page, gid) {
  var url = "/gaea/vertex/find/" + gid;
  window.location.hash = gid

  $.ajax({
    url: url,
    dataType: 'json',
    type: 'GET',
    success: function(result) {
      if (Object.keys(result).length > 0) {
        page.setState({vertex: result, lastMatch: page.state.input})
      }
    }.bind(page),

    error: function(xhr, status, err) {
      console.error(url, status, err.toString());
    }.bind(page)
  });
}

var VertexProperty = React.createClass({
  getInitialState: function() {
    return {};
  },

  render: function() {
    return (
      <div>
        {this.props.name} : {this.props.property}
      </div>
    );
  }
});

var VertexEdge = React.createClass({
  getInitialState: function() {
    return {};
  },

  render: function() {
    return (
      <li onClick={this.props.navigate(this.props.gid)}>{snipPrefix(this.props.gid)}</li>
    )
  }
});

var VertexEdges = React.createClass({
  getInitialState: function() {
    return {};
  },

  render: function() {
    var edges = this;
    var edgeList = this.props.edges.map(function(edge) {
      return <VertexEdge key={edge} gid={edge} navigate={edges.props.navigate} />
    })

    var prefix = this.props.edges[0].split(':')[0]

    return (
      <div>
        <h4>{this.props.label} ({this.props.direction} {prefix})</h4>
        <ul>
          {edgeList}
        </ul>
      </div>
    );
  }
});

var AlternateView = React.createClass({
  getInitialState: function() {
    return {}
  },

  render: function() {
    var page = this;
    return (
      <div onClick={this.props.navigate(this.props.back)}>{this.props.vertex.properties.gid}</div>
    )
  }
});

var VertexView = React.createClass({
  getInitialState: function() {
    return {}
  },

  render: function() {
    var page = this;
    var properties = Object.keys(this.props.vertex.properties).map(function(key) {
      var property = page.props.vertex.properties[key];
      return <VertexProperty key={key} name={key} property={property} />
    });

    var inEdges = Object.keys(this.props.vertex['in']).map(function(key) {
      return <VertexEdges key={key} label={key} navigate={page.props.navigate} edges={page.props.vertex['in'][key]} direction="from"/>
    });

    var outEdges = Object.keys(this.props.vertex['out']).map(function(key) {
      return <VertexEdges key={key} label={key} navigate={page.props.navigate} edges={page.props.vertex['out'][key]} direction="to"/>
    });

    return (
      <div>
        <h2>Showing {this.props.vertex.type} vertex {snipPrefix(this.props.vertex.properties.gid)}</h2>
        <div className="vertex-properties">
          {properties}
        </div>
        <div className="vertex-in-edges">
          <h3>In Edges</h3>
          {inEdges}
        </div>
        <div className="vertex-out-edges">
          <h3>Out Edges</h3>
          {outEdges}
        </div>
      </div>
    )
  }
});

var VertexInput = React.createClass({
  getInitialState: function() {
    var hash = window.location.hash.substr(1);
    exploreVertex(this, hash);

    return {
      input: hash,
      lastMatch: hash,
      back: '',
      vertex: {}
    };
  },

  changeInput: function(event) {
    var gid = event.target.value;
    this.setState({input: gid, back: this.state.input});
    exploreVertex(this, gid);
  },

  render: function() {
    var page = this;
    var navigate = function(gid) {
      return function() {
        page.setState({input: gid, back: page.state.input, lastMatch: gid});
        exploreVertex(page, gid);
      }
    };

    var vertex = <div className="empty-vertex"></div>;
    if (this.state.vertex['properties']) {
      if (this.state.vertex.properties.type === 'Individual') {
        vertex = <VertexView navigate={navigate} vertex={this.state.vertex} back={this.state.back}/>
      } else {
        vertex = <VertexView navigate={navigate} vertex={this.state.vertex} />
      }
    }

    return (
      <div>
        <div className="mdl-textfield mdl-js-textfield mdl-textfield--floating-label">
          <label className="mdl-textfield__label" htmlFor="vertex-gid-input">Enter a vertex GID</label>
          <input id="vertex-gid-input" type="text" name="gid" className="mdl-textfield__input" onChange={(e) => this.changeInput(e)} value={this.state.input} />
        </div>
        {vertex}
      </div>
    );
  }
});

ReactDOM.render(<VertexInput />, document.getElementById('vertex-explore'));
