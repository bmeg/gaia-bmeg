var GraphConsole = React.createClass({
  getInitialState: function() {
    return {
      queries: [{query: "hello", result: "world"}]
    };
  },

  querySubmit: function(e) {
    e.preventDefault()
    e.stopPropagation()

    var query = $('#console-input').val()
    this.state.queries.push({query: query})
    this.setState({queries: this.state.queries})

    $.ajax({
      url: this.props.queryURL,
      dataType: 'json',
      type: 'POST',
      data: JSON.stringify({query: query}),
      success: function(result) {
        var queryTotal = this.state.queries.length
        this.state.queries[queryTotal - 1].result = result.result
        this.setState({queries: this.state.queries});
        $('#console-input').val("")
      }.bind(this),

      error: function(xhr, status, err) {
        console.error(url, status, err.toString());
      }.bind(this)
    });
  },
  
  render: function() {
    var queryLine = function(query) {
      return (
        <div key={query.query+Math.random()} className="mdl-cell mdl-shadow--2dp mdl-cell--12-col-desktop">
          <li className="console-query-line mdl-list__item">
            <span className="mdl-list__item-primary-content">
              <i className="material-icons mdl-list__item-icon mdl-color--white mdl-color-text--teal-100">play_circle_filled</i>
              {query.query}
            </span>
          </li>
          <li className="console-query-line mdl-list__item">
            <span className="mdl-list__item-primary-content">
              {query.result}
            </span>
          </li>
        </div>
      )
    }

    return (
      <div>
        <ul className="mdl-list">
          {this.state.queries.map(queryLine)}
          <li className="console-query-line mdl-list__item">
            <span className="mdl-list__item-primary-content">
              <i className="material-icons mdl-list__item-icon mdl-color--white mdl-color-text--teal-100">play_circle_filled</i>
              <form id="console-submit" onSubmit={(e) => this.querySubmit(e)}>
                <div className="mdl-textfield mdl-js-textfield mdl-textfield--floating-label">
                  <input id="console-input" type="text" name="query" className="mdl-textfield__input" />
                  <label id="console-prompt" className="mdl-textfield__label" htmlFor="console-input">enter query</label>
                </div>
              </form>
            </span>
          </li>
        </ul>
      </div>
    );
  }
});

ReactDOM.render(<GraphConsole queryURL="/gaea/console" />, document.getElementById('graph-console'));
